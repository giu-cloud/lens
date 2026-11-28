import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * La prova del bordo con il codice del motore (28/11): per ogni passaggio e ogni lato, Lens C (Fast con il livello del
 * radar misurato sul riquadro finto) messo in una WindGrid sul solo riquadro finto, e LensWind sopra il modello letto dai
 * punti di mare, come nell'app. Sulle celle della libreria (griglia del radar) scrive la forza (m/s) di:
 *   modello, LensWind di prima (raccordo di 1,5 miglia verso il modello, niente alone), LensWind del motore (con l'alone).
 * Uso: BordoLens <atlante> <libreria> <cartella mappe> <scala km> <rientro km> <uscita>  poi per scena: <GRIB vento> <GRIB profilo> <tempo ms> <u925> <v925>
 * Uscita (little endian): "BRL1", ns, ny, nx della libreria, lat0, lon0, dlat, dlon; per scena: tempo, 4 lati x 3 campi ny*nx float32 (righe da sud).
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
    out.write("BRL1".toByteArray()); wi(scene.size); wi(lib.ny); wi(lib.nx); wd(lib.lat0); wd(lib.lon0); wd(lib.dlat); wd(lib.dlon)
    val nx = atl.nx; val ny = atl.ny
    val S = atl.box[0]; val N = atl.box[1]; val W = atl.box[2]; val E = atl.box[3]
    val dLat = rientro / 111.19; val dLon = rientro / (111.19 * cos(Math.toRadians((S + N) / 2)))
    // righe e colonne dell'atlante (righe da nord) del riquadro finto di ogni lato: r0..r1, c0..c1 inclusi
    val nr = (dLat / atl.cellLat).roundToInt(); val nc = (dLon / atl.cellLon).roundToInt()
    val finti = listOf(intArrayOf(nr, ny - 1, 0, nx - 1), intArrayOf(0, ny - 1 - nr, 0, nx - 1), intArrayOf(0, ny - 1, nc, nx - 1), intArrayOf(0, ny - 1, 0, nx - 1 - nc))
    for ((k, sc) in scene.withIndex()) {
        val (gW, gP, tS, uS, vS) = sc
        val t = tS.toLong()
        val g = WindGrid.from(GribReader.read(Unpack.grib(File(gW).readBytes())))
        if (g == null || t !in g.startMs..g.endMs) { wl(t); repeat(12 * lib.ny * lib.nx) { wf(Float.NaN) }; continue }
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
        wl(t)
        for (fb in finti) {
            val (r0, r1, c0, c1) = fb.toList()
            val bLatN = N - (r0 + 0.5) * atl.cellLat; val bLatS = N - (r1 + 0.5) * atl.cellLat
            val bLonW = W + (c0 + 0.5) * atl.cellLon; val bLonE = W + (c1 + 0.5) * atl.cellLon
            fun dentro(la: Double, lo: Double) = la >= bLatS && la <= bLatN && lo >= bLonW && lo <= bLonE
            val gb = FloatArray(lib.nc) { c -> if (dentro(lib.latCella(c), lib.lonCella(c))) g0[c] else Float.NaN }
            val liv = LensLivello(lib, gb, { la, lo -> if (dentro(la, lo)) forzaL(la, lo) else Double.NaN }, { la, lo -> if (dentro(la, lo)) modello(la, lo) else Double.NaN }, scala)
            // la WindGrid di Lens C sul riquadro finto, righe da sud
            val gnx = c1 - c0 + 1; val gny = r1 - r0 + 1
            val uu = FloatArray(gnx * gny); val vv = FloatArray(gnx * gny)
            for (jj in 0 until gny) for (ii in 0 until gnx) {
                val r = r1 - jj; val c = c0 + ii; val p = r * nx + c
                val lat = N - (r + 0.5) * atl.cellLat; val lon = W + (c + 0.5) * atl.cellLon
                val f = liv.factor(lat, lon).toFloat()
                uu[jj * gnx + ii] = hF.u[p] * f; vv[jj * gnx + ii] = hF.v[p] * f
            }
            val fu = GribField("u", t, t, gnx, gny, bLatS, bLonW, atl.cellLat, atl.cellLon, uu)
            val fv = GribField("v", t, t, gnx, gny, bLatS, bLonW, atl.cellLat, atl.cellLon, vv)
            val grid = WindGrid.from(GribFile(listOf(fu, fv), 4, emptyList()))!!
            val nuovo = LensWind(g, grid)
            val campi = Array(3) { FloatArray(lib.ny * lib.nx) { Float.NaN } }
            for (j in 0 until lib.ny) for (i in 0 until lib.nx) {
                val la = lib.lat0 + j * lib.dlat; val lo = lib.lon0 + i * lib.dlon; val p = j * lib.nx + i
                val wb = g.wind(la, lo, t) ?: continue
                campi[0][p] = (wb.twsKn / KNOTS_PER_MS).toFloat()
                // prima: Lens dentro con il raccordo di 1,5 miglia verso il modello, il modello fuori
                val w = nuovo.weight(la, lo, t)
                val wlv = if (w > 0) grid.wind(la, lo, t) else null
                campi[1][p] = ((if (wlv == null) wb.twsKn else wlv.twsKn * w + wb.twsKn * (1 - w)) / KNOTS_PER_MS).toFloat()
                nuovo.wind(la, lo, t)?.let { campi[2][p] = (it.twsKn / KNOTS_PER_MS).toFloat() }
            }
            for (a in campi) for (x in a) wf(x)
        }
        println("scena ${k + 1}/${scene.size}")
    }
    out.close()
}
