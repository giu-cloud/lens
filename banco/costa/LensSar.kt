import app.regata.core.*
import java.io.*
import java.nio.*
import kotlin.math.*

/** Lens sulle scene del radar: per ogni scena il vento del modello (ECMWF nel prodotto SAR) passa per l'atlante, e il
 *  risultato si riporta sulla griglia del radar. Ingresso e uscita: file binari scritti e letti da costa_analisi.py. */
fun main(args: Array<String>) {
    val atlas = LensAtlas.read(File(args[0]).nameWithoutExtension, File(args[0]).readBytes())
    val inp = DataInputStream(BufferedInputStream(FileInputStream(args[1])))
    fun i32() = Integer.reverseBytes(inp.readInt())
    fun f64() = java.lang.Double.longBitsToDouble(java.lang.Long.reverseBytes(inp.readLong()))
    fun i64() = java.lang.Long.reverseBytes(inp.readLong())
    fun f32() = java.lang.Float.intBitsToFloat(Integer.reverseBytes(inp.readInt()))
    val magic = ByteArray(4).also { inp.readFully(it) }
    require(String(magic) == "SAR1")
    val ny = i32(); val nx = i32(); val lat0 = f64(); val dlat = f64(); val lon0 = f64(); val dlon = f64(); val ns = i32()
    val n = (System.getProperty("lens.n") ?: "${LensAtlas.N_MIXED}").toDouble()
    val out = DataOutputStream(BufferedOutputStream(FileOutputStream(args[2])))
    fun w32(x: Int) = out.writeInt(Integer.reverseBytes(x))
    fun wf(x: Float) = out.writeInt(Integer.reverseBytes(java.lang.Float.floatToIntBits(x)))
    out.write("LSO1".toByteArray()); w32(ny); w32(nx); w32(ns)
    val box = atlas.box
    val bnx = atlas.nx; val bny = atlas.ny
    val cLat = atlas.cellLat; val cLon = atlas.cellLon
    var done = 0
    for (s in 0 until ns) {
        val t = i64()
        val mu = FloatArray(ny * nx) { f32() }; val mv = FloatArray(ny * nx) { f32() }
        // il vento del modello riempito dove manca (terra, bordi) con la media della scena
        var su = 0.0; var sv = 0.0; var c = 0
        for (k in mu.indices) if (mu[k].isFinite() && mv[k].isFinite()) { su += mu[k]; sv += mv[k]; c++ }
        if (c < 20) { for (k in 0 until ny * nx * 2) wf(Float.NaN); continue }
        val mU = su / c; val mV = sv / c
        val fu = FloatArray(ny * nx) { if (mu[it].isFinite()) mu[it] else mU.toFloat() }
        val fv = FloatArray(ny * nx) { if (mv[it].isFinite()) mv[it] else mV.toFloat() }
        val base = object : WindProvider {
            override val startMs = t - 3_600_000L
            override val endMs = t + 3_600_000L
            override fun wind(lat: Double, lon: Double, timeMs: Long): Wind {
                val y = ((lat - lat0) / dlat).coerceIn(0.0, ny - 1.0); val x = ((lon - lon0) / dlon).coerceIn(0.0, nx - 1.0)
                val y0 = floor(y).toInt().coerceAtMost(ny - 2).coerceAtLeast(0); val x0 = floor(x).toInt().coerceAtMost(nx - 2).coerceAtLeast(0)
                val ty = y - y0; val tx = x - x0
                fun bl(a: FloatArray) = (a[y0 * nx + x0] * (1 - tx) + a[y0 * nx + x0 + 1] * tx) * (1 - ty) + (a[(y0 + 1) * nx + x0] * (1 - tx) + a[(y0 + 1) * nx + x0 + 1] * tx) * ty
                val u = bl(fu); val v = bl(fv)
                val sp = hypot(u, v)
                return Wind(Angles.norm360(Math.toDegrees(atan2(-u, -v))), sp * KNOTS_PER_MS)
            }
        }
        val hr = atlas.hour(base, null, t, n)
        for (comp in 0..1) {
            val a = if (hr == null) null else if (comp == 0) hr.u else hr.v
            for (j in 0 until ny) for (i in 0 until nx) {
                if (a == null) { wf(Float.NaN); continue }
                val lat = lat0 + j * dlat; val lon = lon0 + i * dlon
                // griglia di Lens: righe da nord, centri delle celle
                val fy = (box[1] - lat) / cLat - 0.5; val fx = (lon - box[2]) / cLon - 0.5
                if (fy < -0.5 || fy > bny - 0.5 || fx < -0.5 || fx > bnx - 0.5) { wf(Float.NaN); continue }
                val yy = fy.coerceIn(0.0, bny - 1.0); val xx = fx.coerceIn(0.0, bnx - 1.0)
                val y0 = floor(yy).toInt().coerceAtMost(bny - 2); val x0 = floor(xx).toInt().coerceAtMost(bnx - 2)
                val ty = yy - y0; val tx = xx - x0
                val v = (a[y0 * bnx + x0] * (1 - tx) + a[y0 * bnx + x0 + 1] * tx) * (1 - ty) + (a[(y0 + 1) * bnx + x0] * (1 - tx) + a[(y0 + 1) * bnx + x0 + 1] * tx) * ty
                wf(v.toFloat())
            }
        }
        done++
    }
    out.close()
    System.err.println("scene con Lens: $done su $ns")
}
