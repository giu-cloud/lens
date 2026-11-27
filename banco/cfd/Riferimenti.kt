import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * I riferimenti della prova CFD (27/11, lens/cfd_piano_26_11.md): per i passaggi del radar di un'area, sull'area grande,
 * il vento come lo darebbe l'app con la stabilità vera della scena (il profilo ECMWF del GRIB "stab"):
 *  - M: il modello letto dai soli punti di mare (come Shift); M10: il modello lisciato come lo riceve Lens;
 *  - Fast: l'atlante (LensAtlas.hour) con l'aria della scena; Detailed: l'ora risolta (LensHourly.fields) con la stessa aria;
 *  - il livello del radar (LensLivello, scala dell'area) misurato sull'uno e sull'altro;
 *  - il fattore di B (relativoALens sul Fast), come riferimento di quello che il radar sa.
 * Con -Dshift.lens.relief=k la manopola del rilievo del motore (1 = A, 0 = la 0.54; la rotazione resta quella piena);
 * con -Dlens.soloFast=true niente Detailed (copia di Fast).
 * Uso: Riferimenti <atlante> <libreria> <cartella mappe> <scala km> <uscita>  poi per scena: <GRIB vento> <GRIB profilo> <tempo ms> <u925> <v925>
 * Uscita (little endian): "RIF1", ns, nx, ny, box (4 double); per scena: tempo (long), Fr, N, 12 campi nx*ny float32 a righe da nord:
 *   M, M10, uFast, vFast, uDet, vDet, livFast, livDet, fattore B, mare, Fr del Detailed sotto le cime (costante), secondi del Detailed.
 */
/**
 * L'aria della scena dal GRIB "stab" di ECMWF (t a 1000, 925 e 850 hPa, u e v a 925 e 850, pressione al mare, t a 2 m):
 * come LensAir.fromGrid (media dei 9 punti dell'area), ma le quote dei livelli vengono dall'ipsometrica perché il GRIB non
 * ha il geopotenziale, e il vento a 1000 hPa è 1,2 volte quello a 10 m del modello.
 */
fun aria(f: File, g: WindGrid, atl: LensAtlas, t: Long): LensAirHour? {
    val file = runCatching { GribReader.read(Unpack.grib(f.readBytes())) }.getOrNull() ?: return null
    fun sg(n: String) = ScalarGrid.from(file, n)
    val T = listOf("t1000", "t925", "t850").map { sg(it) ?: return null }
    val U = listOf("u925", "u850").map { sg(it) ?: return null }; val V = listOf("v925", "v850").map { sg(it) ?: return null }
    val P = sg("prmsl") ?: return null; val T2 = sg("t2m") ?: return null
    val lv = intArrayOf(1000, 925, 850)
    val Z = DoubleArray(3); val Tc = DoubleArray(3); val Us = DoubleArray(3); val Ux = DoubleArray(3); val Vx = DoubleArray(3); var n = 0
    for ((la, lo) in LensAir.points(atl.domain)) {
        val tk = T.map { it.value(la, lo, t) ?: Double.NaN }; val p0 = P.value(la, lo, t) ?: continue; val t2 = T2.value(la, lo, t) ?: continue
        if (tk.any { !it.isFinite() }) continue
        val z1000 = 29.27 * 0.5 * (t2 + tk[0]) * ln(p0 / 100000.0)
        val z925 = z1000 + 29.27 * 0.5 * (tk[0] + tk[1]) * ln(1000.0 / 925.0)
        val z850 = z925 + 29.27 * 0.5 * (tk[1] + tk[2]) * ln(925.0 / 850.0)
        val w = g.wind(la, lo, t)
        val u10 = if (w != null) -w.twsKn / KNOTS_PER_MS * sin(Math.toRadians(w.twdDeg)) * 1.2 else Double.NaN
        val v10 = if (w != null) -w.twsKn / KNOTS_PER_MS * cos(Math.toRadians(w.twdDeg)) * 1.2 else Double.NaN
        val uu = doubleArrayOf(u10, U[0].value(la, lo, t) ?: Double.NaN, U[1].value(la, lo, t) ?: Double.NaN)
        val vv = doubleArrayOf(v10, V[0].value(la, lo, t) ?: Double.NaN, V[1].value(la, lo, t) ?: Double.NaN)
        if (uu.any { !it.isFinite() } || vv.any { !it.isFinite() }) continue
        val zz = doubleArrayOf(z1000, z925, z850)
        for (k in 0 until 3) { Z[k] += zz[k]; Tc[k] += tk[k] - 273.15; Ux[k] += uu[k]; Vx[k] += vv[k]; Us[k] += hypot(uu[k], vv[k]) }
        n++
    }
    if (n == 0) return null
    for (k in 0 until 3) { Z[k] /= n; Tc[k] /= n; Ux[k] /= n; Vx[k] /= n; Us[k] /= n }
    return LensAir.hourFrom(lv, Z, Tc, Us, Ux, Vx, BooleanArray(3) { true }, atl.hRil)
}

fun main(args: Array<String>) {
    val atl = LensAtlas.read(File(args[0]).nameWithoutExtension.removePrefix("atlante_"), File(args[0]).readBytes())
    val lib = LensRadarLibreria.leggi(File(args[1]))
    val maps = args[2]; val scala = args[3].toDouble()
    val coast = CoastIndex(MapPack.read(File("$maps/europa.rgm").inputStream()), MapPack.read(File("$maps/world.rgm").inputStream()), listOf(MapPack.read(File("$maps/atlantico.rgm").inputStream())))
    WindGrid.landTest = { la, lo -> coast.onLand(la, lo) }
    val scene = args.drop(5).chunked(5)
    val out = DataOutputStream(BufferedOutputStream(FileOutputStream(args[4])))
    fun wi(x: Int) = out.writeInt(Integer.reverseBytes(x))
    fun wl(x: Long) = out.writeLong(java.lang.Long.reverseBytes(x))
    fun wd(x: Double) = out.writeLong(java.lang.Long.reverseBytes(java.lang.Double.doubleToLongBits(x)))
    fun wf(x: Float) = out.writeInt(Integer.reverseBytes(java.lang.Float.floatToIntBits(x)))
    out.write("RIF1".toByteArray()); wi(scene.size); wi(atl.nx); wi(atl.ny); for (b in atl.box) wd(b)
    val nx = atl.nx; val ny = atl.ny
    val sea = atl.seaInner()
    for ((k, sc) in scene.withIndex()) {
        val (gW, gP, tS, uS, vS) = sc
        val t = tS.toLong()
        val g = WindGrid.from(GribReader.read(Unpack.grib(File(gW).readBytes())))
        val fld = Array(12) { FloatArray(nx * ny) { Float.NaN } }
        if (g == null || t !in g.startMs..g.endMs) {
            System.err.println("scena $k: il GRIB del vento non copre l'ora"); wl(t); wd(Double.NaN); wd(Double.NaN); for (a in fld) for (x in a) wf(x); continue
        }
        val air = aria(File(gP), g, atl, t)
        // Fast e Detailed con l'aria della scena
        val hF = atl.hour(g, null, t, air?.n ?: LensAtlas.N_MIXED, air)!!
        val t0 = System.nanoTime()
        val uD = FloatArray(nx * ny); val vD = FloatArray(nx * ny)
        if (System.getProperty("lens.soloFast") == "true") {          // Detailed = Fast (per le prove della manopola del rilievo)
            hF.u.copyInto(uD); hF.v.copyInto(vD)
        } else {
            val det = LensHourly.fields(atl, g, null, listOf(t), air = { air }, n = air?.n ?: LensAtlas.N_MIXED, threads = 2)
            val dU = det.first { it.param == "u" }.values; val dV = det.first { it.param == "v" }.values   // righe da sud
            for (r in 0 until ny) for (c in 0 until nx) { uD[r * nx + c] = dU[(ny - 1 - r) * nx + c]; vD[r * nx + c] = dV[(ny - 1 - r) * nx + c] }
        }
        val sec = (System.nanoTime() - t0) / 1e9
        // il modello sulle celle della libreria e sull'atlante (righe da sud, come Lens2Use.modelOnAtlas)
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
        val modello = atl.speedSampler(bU, bV, rowsFromSouth = true)
        val livF = LensLivello(lib, g0, atl.speedSampler(hF.u, hF.v), modello, scala)
        val livD = LensLivello(lib, g0, atl.speedSampler(uD, vD), modello, scala)
        val campoB = LensRadarCampo(lib, LensRadar.relativoALens(lib, g0, FloatArray(lib.nc) { hypot(mu[it], mv[it]) }, atl.speedSampler(hF.u, hF.v)))
        val m10 = atl.model10(g, null, t)!!
        for (r in 0 until ny) {
            val lat = atl.box[1] - (r + 0.5) * atl.cellLat
            for (c in 0 until nx) {
                val lon = atl.box[2] + (c + 0.5) * atl.cellLon; val p = r * nx + c
                vento(lat, lon)?.let { fld[0][p] = hypot(it.first, it.second) }
                fld[1][p] = m10.s[(r + atl.margin) * atl.NX + c + atl.margin].toFloat()
                fld[2][p] = hF.u[p]; fld[3][p] = hF.v[p]; fld[4][p] = uD[p]; fld[5][p] = vD[p]
                fld[6][p] = livF.factor(lat, lon).toFloat(); fld[7][p] = livD.factor(lat, lon).toFloat()
                fld[8][p] = campoB.fattore(lat, lon); fld[9][p] = if (sea[p]) 1f else 0f
                fld[10][p] = (air?.fr ?: Double.NaN).toFloat(); fld[11][p] = sec.toFloat()
            }
        }
        wl(t); wd(air?.fr ?: Double.NaN); wd(air?.n ?: Double.NaN)
        for (a in fld) for (x in a) wf(x)
        println("scena ${k + 1}/${scene.size} ${java.time.Instant.ofEpochMilli(t)}: Fr ${"%.2f".format(air?.fr ?: Double.NaN)}, Detailed ${"%.1f".format(sec)} s")
    }
    out.close()
}
