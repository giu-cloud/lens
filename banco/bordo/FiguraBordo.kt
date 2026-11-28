import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * Figura del bordo (28/11): Lens C sull'area grande intera, e attorno, su una griglia più larga di [margine] gradi a passo
 * [passo]: la forza e la direzione (u, v in m/s) del modello, di LensWind di prima (raccordo di 1,5 miglia, niente alone) e
 * di LensWind del motore (con l'alone). Uso: FiguraBordo <atlante> <libreria> <mappe> <scala km> <margine> <passo> <uscita> <GRIB> <profilo> <tempo ms>
 * Uscita: "FBD1", ny, nx, lat0, lon0, passo, box (4 double), poi 6 campi ny*nx float32 (righe da sud): uM, vM, uPrima, vPrima, uAlone, vAlone.
 */
fun main(a: Array<String>) {
    val atl = LensAtlas.read(File(a[0]).nameWithoutExtension.removePrefix("atlante_"), File(a[0]).readBytes())
    val lib = LensRadarLibreria.leggi(File(a[1])); val maps = a[2]; val scala = a[3].toDouble(); val mg = a[4].toDouble(); val ps = a[5].toDouble()
    val coast = CoastIndex(MapPack.read(File("$maps/europa.rgm").inputStream()), MapPack.read(File("$maps/world.rgm").inputStream()), listOf(MapPack.read(File("$maps/atlantico.rgm").inputStream())))
    WindGrid.landTest = { la, lo -> coast.onLand(la, lo) }
    val t = a[9].toLong(); val g = WindGrid.from(GribReader.read(Unpack.grib(File(a[7]).readBytes())))!!
    val air = aria(File(a[8]), g, atl, t)
    val hF = atl.hour(g, null, t, air?.n ?: LensAtlas.N_MIXED, air)!!
    fun vento(la: Double, lo: Double): Pair<Float, Float>? { val w = g.wind(la, lo, t) ?: return null; val ms = w.twsKn / KNOTS_PER_MS; val r = Math.toRadians(w.twdDeg); return (-ms * sin(r)).toFloat() to (-ms * cos(r)).toFloat() }
    val nx = atl.nx; val ny = atl.ny
    val mu = FloatArray(lib.nc) { Float.NaN }; val mv = FloatArray(lib.nc) { Float.NaN }
    for (c in 0 until lib.nc) vento(lib.latCella(c), lib.lonCella(c))?.let { mu[c] = it.first; mv[c] = it.second }
    val bU = FloatArray(nx * ny) { Float.NaN }; val bV = FloatArray(nx * ny) { Float.NaN }
    for (j in 0 until ny) for (i in 0 until nx) vento(atl.box[0] + (j + 0.5) * atl.cellLat, atl.box[2] + (i + 0.5) * atl.cellLon)?.let { bU[j * nx + i] = it.first; bV[j * nx + i] = it.second }
    val ok = (0 until lib.nc).filter { !mu[it].isNaN() }
    val s = LensRadar.situazione(lib, mu, mv, ok.sumOf { mu[it].toDouble() } / ok.size, ok.sumOf { mv[it].toDouble() } / ok.size, ((t / 3_600_000L) % 24L).toInt())
    val g0 = lib.fattore(s, FloatArray(lib.nc) { LensRadar.provenienza(mu[it], mv[it]) })
    val liv = LensLivello(lib, g0, atl.speedSampler(hF.u, hF.v), atl.speedSampler(bU, bV, rowsFromSouth = true), scala)
    val uu = FloatArray(nx * ny); val vv = FloatArray(nx * ny)
    for (jj in 0 until ny) for (i in 0 until nx) {
        val r = ny - 1 - jj; val p = r * nx + i
        val f = liv.factor(atl.box[1] - (r + 0.5) * atl.cellLat, atl.box[2] + (i + 0.5) * atl.cellLon).toFloat()
        uu[jj * nx + i] = hF.u[p] * f; vv[jj * nx + i] = hF.v[p] * f
    }
    val grid = WindGrid.from(GribFile(listOf(GribField("u", t, t, nx, ny, atl.box[0] + 0.5 * atl.cellLat, atl.box[2] + 0.5 * atl.cellLon, atl.cellLat, atl.cellLon, uu),
        GribField("v", t, t, nx, ny, atl.box[0] + 0.5 * atl.cellLat, atl.box[2] + 0.5 * atl.cellLon, atl.cellLat, atl.cellLon, vv)), 4, emptyList()))!!
    val lw = LensWind(g, grid)
    val la0 = atl.box[0] - mg; val lo0 = atl.box[2] - mg
    val NY = ((atl.box[1] - atl.box[0] + 2 * mg) / ps).toInt() + 1; val NX = ((atl.box[3] - atl.box[2] + 2 * mg) / ps).toInt() + 1
    val F = Array(6) { FloatArray(NY * NX) { Float.NaN } }
    fun uv(w: Wind?): Pair<Float, Float>? { w ?: return null; val ms = w.twsKn / KNOTS_PER_MS; val r = Math.toRadians(w.twdDeg); return (-ms * sin(r)).toFloat() to (-ms * cos(r)).toFloat() }
    for (j in 0 until NY) for (i in 0 until NX) {
        val la = la0 + j * ps; val lo = lo0 + i * ps; val p = j * NX + i
        val wb = g.wind(la, lo, t) ?: continue
        uv(wb)?.let { F[0][p] = it.first; F[1][p] = it.second }
        val w = lw.weight(la, lo, t); val wl = if (w > 0) grid.wind(la, lo, t) else null
        val prima = if (wl == null) wb else {
            val al = Math.toRadians(wl.twdDeg); val bl = Math.toRadians(wb.twdDeg)
            val u = -(wl.twsKn * sin(al)) * w - (wb.twsKn * sin(bl)) * (1 - w); val v = -(wl.twsKn * cos(al)) * w - (wb.twsKn * cos(bl)) * (1 - w)
            Wind(Angles.norm360(Math.toDegrees(atan2(-u, -v))), wl.twsKn * w + wb.twsKn * (1 - w))
        }
        uv(prima)?.let { F[2][p] = it.first; F[3][p] = it.second }
        uv(lw.wind(la, lo, t))?.let { F[4][p] = it.first; F[5][p] = it.second }
    }
    val out = DataOutputStream(BufferedOutputStream(FileOutputStream(a[6])))
    fun wi(x: Int) = out.writeInt(Integer.reverseBytes(x))
    fun wd(x: Double) = out.writeLong(java.lang.Long.reverseBytes(java.lang.Double.doubleToLongBits(x)))
    out.write("FBD1".toByteArray()); wi(NY); wi(NX); wd(la0); wd(lo0); wd(ps); for (b in atl.box) wd(b)
    for (f in F) for (x in f) out.writeInt(Integer.reverseBytes(java.lang.Float.floatToIntBits(x)))
    out.close()
}
