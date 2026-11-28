import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * Il bordo dei riquadri di Lens (28/11, lens/bordo_28_11.md). Per i passaggi del radar di un'area grande, sull'atlante:
 * il modello come lo riceve Lens (M10), la forza di Lens Fast (manopola del motore) e il livello del radar (LensLivello)
 * misurato su tutta l'area e, per ciascuno dei quattro lati, su un riquadro "finto" più piccolo: il lato spostato verso
 * l'interno di [km] chilometri, con il radar e Lens fuori dal riquadro finto tolti (come se l'area finisse lì).
 * Uso: Bordo <atlante> <libreria> <cartella mappe> <scala km> <rientro km> <uscita>  poi per scena: <GRIB vento> <GRIB profilo> <tempo ms> <u925> <v925>
 * Uscita (little endian): "BRD1", ns, nx, ny, box (4 double); per scena: tempo (long), 8 campi nx*ny float32 a righe da nord:
 *   M10, forza Fast, mare, livello pieno, livello con il lato nord rientrato, sud, ovest, est.
 */
fun main(args: Array<String>) {
    val atl = LensAtlas.read(File(args[0]).nameWithoutExtension.removePrefix("atlante_"), File(args[0]).readBytes())
    val lib = LensRadarLibreria.leggi(File(args[1]))
    val maps = args[2]; val scala = args[3].toDouble(); val rientro = args[4].toDouble()
    val coast = CoastIndex(MapPack.read(File("$maps/europa.rgm").inputStream()), MapPack.read(File("$maps/world.rgm").inputStream()), listOf(MapPack.read(File("$maps/atlantico.rgm").inputStream())))
    WindGrid.landTest = { la, lo -> coast.onLand(la, lo) }
    val scene = args.drop(6).chunked(5)
    val out = DataOutputStream(BufferedOutputStream(FileOutputStream(args[5])))
    fun wi(x: Int) = out.writeInt(Integer.reverseBytes(x))
    fun wl(x: Long) = out.writeLong(java.lang.Long.reverseBytes(x))
    fun wd(x: Double) = out.writeLong(java.lang.Long.reverseBytes(java.lang.Double.doubleToLongBits(x)))
    fun wf(x: Float) = out.writeInt(Integer.reverseBytes(java.lang.Float.floatToIntBits(x)))
    out.write("BRD1".toByteArray()); wi(scene.size); wi(atl.nx); wi(atl.ny); for (b in atl.box) wd(b)
    val nx = atl.nx; val ny = atl.ny
    val sea = atl.seaInner()
    val S = atl.box[0]; val N = atl.box[1]; val W = atl.box[2]; val E = atl.box[3]
    val dLat = rientro / 111.19; val dLon = rientro / (111.19 * cos(Math.toRadians((S + N) / 2)))
    // i quattro riquadri finti: nord, sud, ovest, est rientrati
    val finti = listOf(doubleArrayOf(S, N - dLat, W, E), doubleArrayOf(S + dLat, N, W, E), doubleArrayOf(S, N, W + dLon, E), doubleArrayOf(S, N, W, E - dLon))
    fun dentro(b: DoubleArray, la: Double, lo: Double) = la >= b[0] && la <= b[1] && lo >= b[2] && lo <= b[3]
    for ((k, sc) in scene.withIndex()) {
        val (gW, gP, tS, uS, vS) = sc
        val t = tS.toLong()
        val g = WindGrid.from(GribReader.read(Unpack.grib(File(gW).readBytes())))
        val fld = Array(8) { FloatArray(nx * ny) { Float.NaN } }
        if (g == null || t !in g.startMs..g.endMs) { wl(t); for (a in fld) for (x in a) wf(x); continue }
        val air = aria(File(gP), g, atl, t)
        val hF = atl.hour(g, null, t, air?.n ?: LensAtlas.N_MIXED, air)!!
        fun vento(la: Double, lo: Double): Pair<Float, Float>? {
            val w = g.wind(la, lo, t) ?: return null
            val ms = w.twsKn / KNOTS_PER_MS; val a = Math.toRadians(w.twdDeg)
            return (-ms * sin(a)).toFloat() to (-ms * cos(a)).toFloat()
        }
        val mu = FloatArray(lib.nc) { Float.NaN }; val mv = FloatArray(lib.nc) { Float.NaN }
        for (c in 0 until lib.nc) vento(lib.latCella(c), lib.lonCella(c))?.let { mu[c] = it.first; mv[c] = it.second }
        val bU = FloatArray(nx * ny) { Float.NaN }; val bV = FloatArray(nx * ny) { Float.NaN }
        for (j in 0 until ny) for (i in 0 until nx) vento(atl.box[0] + (j + 0.5) * atl.cellLat, atl.box[2] + (i + 0.5) * atl.cellLon)?.let { bU[j * nx + i] = it.first; bV[j * nx + i] = it.second }
        var u925 = uS.toDoubleOrNull() ?: Double.NaN; var v925 = vS.toDoubleOrNull() ?: Double.NaN
        if (u925.isNaN() || v925.isNaN()) { val ok = (0 until lib.nc).filter { !mu[it].isNaN() }; u925 = ok.sumOf { mu[it].toDouble() } / ok.size; v925 = ok.sumOf { mv[it].toDouble() } / ok.size }
        val s = LensRadar.situazione(lib, mu, mv, u925, v925, ((t / 3_600_000L) % 24L).toInt())
        val g0 = lib.fattore(s, FloatArray(lib.nc) { LensRadar.provenienza(mu[it], mv[it]) })
        val forzaL = atl.speedSampler(hF.u, hF.v)
        val modello = atl.speedSampler(bU, bV, rowsFromSouth = true)
        val livelli = mutableListOf(LensLivello(lib, g0, forzaL, modello, scala))
        for (b in finti) {
            val gb = FloatArray(lib.nc) { c -> if (dentro(b, lib.latCella(c), lib.lonCella(c))) g0[c] else Float.NaN }
            livelli += LensLivello(lib, gb, { la, lo -> if (dentro(b, la, lo)) forzaL(la, lo) else Double.NaN }, { la, lo -> if (dentro(b, la, lo)) modello(la, lo) else Double.NaN }, scala)
        }
        val m10 = atl.model10(g, null, t)!!
        for (r in 0 until ny) {
            val lat = atl.box[1] - (r + 0.5) * atl.cellLat
            for (c in 0 until nx) {
                val lon = atl.box[2] + (c + 0.5) * atl.cellLon; val p = r * nx + c
                fld[0][p] = m10.s[(r + atl.margin) * atl.NX + c + atl.margin].toFloat()
                fld[1][p] = hypot(hF.u[p], hF.v[p]); fld[2][p] = if (sea[p]) 1f else 0f
                for (q in 0 until 5) fld[3 + q][p] = livelli[q].factor(lat, lon).toFloat()
            }
        }
        wl(t); for (a in fld) for (x in a) wf(x)
        println("scena ${k + 1}/${scene.size}")
    }
    out.close()
}
