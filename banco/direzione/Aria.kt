import app.regata.core.*
import java.io.*
import kotlin.math.*

/** L'aria della scena dal GRIB "stab" di ECMWF: copia di Riferimenti.kt (banco/cfd). */
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

