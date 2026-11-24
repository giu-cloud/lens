import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * Il vento di ECMWF come lo legge Shift (interpolazione bilineare dei punti della griglia, terra compresa) e come lo
 * leggerebbe usando solo i punti di mare, sulle celle del radar; e Lens su tutti e due.
 * Ingresso (scritto da costa_ecmwf.py): "SAR2", ny, nx, lat0, dlat, lon0, dlon, ns, poi per scena: tempo (ms) e il nome
 * del GRIB. Uscita "SHO1": per scena 8 campi ny*nx float32 (u,v Shift; u,v solo mare; u,v Lens su Shift; u,v Lens su mare).
 */
val MASK: String = System.getProperty("sea.mask") ?: "point"
val FRAC: Double = System.getProperty("sea.frac")?.toDoubleOrNull() ?: 0.5
val LAND_W: Double = System.getProperty("sea.landw")?.toDoubleOrNull() ?: 0.0

fun main(args: Array<String>) {
    val assets = "/home/claude/Routing/app/src/main/assets"
    val coast = CoastIndex(MapPack.read(File("$assets/maps/europa.rgm").inputStream()), MapPack.read(File("$assets/maps/world.rgm").inputStream()), listOf(MapPack.read(File("$assets/maps/atlantico.rgm").inputStream())))
    val atlas = if (args[0] == "-") null else LensAtlas.read(File(args[0]).nameWithoutExtension, File(args[0]).readBytes())
    val inp = DataInputStream(BufferedInputStream(FileInputStream(args[1])))
    val gribDir = args[3]
    fun i32() = Integer.reverseBytes(inp.readInt())
    fun f64() = java.lang.Double.longBitsToDouble(java.lang.Long.reverseBytes(inp.readLong()))
    fun i64() = java.lang.Long.reverseBytes(inp.readLong())
    require(String(ByteArray(4).also { inp.readFully(it) }) == "SAR2")
    val ny = i32(); val nx = i32(); val lat0 = f64(); val dlat = f64(); val lon0 = f64(); val dlon = f64(); val ns = i32()
    val out = DataOutputStream(BufferedOutputStream(FileOutputStream(args[2])))
    fun wf(x: Float) = out.writeInt(Integer.reverseBytes(java.lang.Float.floatToIntBits(x)))
    out.write("SHO1".toByteArray()); out.writeInt(Integer.reverseBytes(ny)); out.writeInt(Integer.reverseBytes(nx)); out.writeInt(Integer.reverseBytes(ns))
    val cache = HashMap<String, WindGrid?>()
    val landCache = HashMap<String, BooleanArray>()
    var done = 0
    for (s in 0 until ns) {
        val t = i64()
        val nb = i32(); val name = String(ByteArray(nb).also { inp.readFully(it) })
        val g = cache.getOrPut(name) { runCatching { WindGrid.from(GribReader.read(Unpack.grib(File(gribDir, name).readBytes()))) }.getOrNull() }
        if (g == null || t < g.startMs || t > g.endMs) { repeat(8 * ny * nx) { wf(Float.NaN) }; continue }
        // i punti della griglia del modello che stanno a terra
        val land = landCache.getOrPut("${g.ni}x${g.nj}:${g.lat0}:${g.lon0}") { BooleanArray(g.ni * g.nj) { k -> val j = k / g.ni; val i = k % g.ni
            val la = g.lat0 + j * g.dLat; val lo = g.lon0 + i * g.dLon
            if (MASK == "point") coast.onLand(la, lo) else {
                // la cella del modello attorno al punto: terra se la sua parte di terra supera la soglia (come la maschera del modello)
                var n = 0; for (a in -2..2) for (b in -2..2) if (coast.onLand(la + a * g.dLat / 5, lo + b * g.dLon / 5)) n++
                n / 25.0 > FRAC
            } } }
        fun uvAt(lat: Double, lon: Double): DoubleArray? { val w = g.wind(lat, lon, t) ?: return null; val sp = w.twsKn / KNOTS_PER_MS; val a = Math.toRadians(w.twdDeg); return doubleArrayOf(-sp * sin(a), -sp * cos(a)) }
        val nodeU = HashMap<Int, DoubleArray?>()
        fun node(j: Int, i: Int): DoubleArray? { if (j !in 0 until g.nj || i !in 0 until g.ni) return null; return nodeU.getOrPut(j * g.ni + i) { uvAt(g.lat0 + j * g.dLat, g.lon0 + i * g.dLon) } }
        /** Solo i punti di mare: i quattro vicini che sono in mare, pesati come nella bilineare; se nessuno, il più vicino in mare entro 2 celle. */
        fun seaUv(lat: Double, lon: Double): DoubleArray? {
            val fy = (lat - g.lat0) / g.dLat; val fx = (lon - g.lon0) / g.dLon
            val j0 = floor(fy).toInt(); val i0 = floor(fx).toInt(); val ty = fy - j0; val tx = fx - i0
            var su = 0.0; var sv = 0.0; var sw = 0.0
            for ((dj, di, w) in listOf(Triple(0, 0, (1 - ty) * (1 - tx)), Triple(0, 1, (1 - ty) * tx), Triple(1, 0, ty * (1 - tx)), Triple(1, 1, ty * tx))) {
                val j = j0 + dj; val i = i0 + di
                if (j !in 0 until g.nj || i !in 0 until g.ni) continue
                val ww = if (land[j * g.ni + i]) w * LAND_W else w
                if (ww <= 0.0) continue
                val v = node(j, i) ?: continue
                su += v[0] * ww; sv += v[1] * ww; sw += ww
            }
            if (sw > 1e-6) return doubleArrayOf(su / sw, sv / sw)
            var best: DoubleArray? = null; var bd = Double.MAX_VALUE
            for (j in j0 - 2..j0 + 3) for (i in i0 - 2..i0 + 3) {
                if (j !in 0 until g.nj || i !in 0 until g.ni || land[j * g.ni + i]) continue
                val d = (j - fy) * (j - fy) + (i - fx) * (i - fx)
                if (d < bd) { bd = d; best = node(j, i) }
            }
            return best ?: uvAt(lat, lon)
        }
        val fields = Array(8) { FloatArray(ny * nx) { Float.NaN } }
        for (j in 0 until ny) for (i in 0 until nx) {
            val lat = lat0 + j * dlat; val lon = lon0 + i * dlon
            uvAt(lat, lon)?.let { fields[0][j * nx + i] = it[0].toFloat(); fields[1][j * nx + i] = it[1].toFloat() }
            seaUv(lat, lon)?.let { fields[2][j * nx + i] = it[0].toFloat(); fields[3][j * nx + i] = it[1].toFloat() }
        }
        if (atlas != null) {
            val seaProvider = object : WindProvider {
                override val startMs = g.startMs; override val endMs = g.endMs
                override fun wind(lat: Double, lon: Double, timeMs: Long): Wind? { val v = seaUv(lat, lon) ?: return null; val sp = hypot(v[0], v[1]); return Wind(Angles.norm360(Math.toDegrees(atan2(-v[0], -v[1]))), sp * KNOTS_PER_MS) }
            }
            for ((k, prov) in listOf(4 to (g as WindProvider), 6 to (seaProvider as WindProvider))) {
                val hr = atlas.hour(prov, null, t) ?: continue
                val box = atlas.box; val bnx = atlas.nx; val bny = atlas.ny
                for (j in 0 until ny) for (i in 0 until nx) {
                    val lat = lat0 + j * dlat; val lon = lon0 + i * dlon
                    val fy = (box[1] - lat) / atlas.cellLat - 0.5; val fx = (lon - box[2]) / atlas.cellLon - 0.5
                    if (fy < -0.5 || fy > bny - 0.5 || fx < -0.5 || fx > bnx - 0.5) continue
                    val yy = fy.coerceIn(0.0, bny - 1.0); val xx = fx.coerceIn(0.0, bnx - 1.0)
                    val y0 = floor(yy).toInt().coerceAtMost(bny - 2); val x0 = floor(xx).toInt().coerceAtMost(bnx - 2)
                    val ty = yy - y0; val tx = xx - x0
                    fun bl(a: FloatArray) = (a[y0 * bnx + x0] * (1 - tx) + a[y0 * bnx + x0 + 1] * tx) * (1 - ty) + (a[(y0 + 1) * bnx + x0] * (1 - tx) + a[(y0 + 1) * bnx + x0 + 1] * tx) * ty
                    fields[k][j * nx + i] = bl(hr.u).toFloat(); fields[k + 1][j * nx + i] = bl(hr.v).toFloat()
                }
            }
        }
        for (f in fields) for (x in f) wf(x)
        done++
    }
    out.close()
    System.err.println("scene con ECMWF: $done su $ns")
}
