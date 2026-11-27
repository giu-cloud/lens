import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * Come Schermate.kt, per qualsiasi area (27/11 sera, revisione): atlante grande, atlante fine, libreria, mappe, uscita,
 * scala del livello A e B (km), poi per scena <GRIB> <tempo ms> <u925> <v925>. Campi come Schermate.kt, con il fattore del
 * livello alla scala A nel campo 6 e alla scala B nel campo 7. La manopola del rilievo è quella del motore (o -Dshift.lens.relief).
 */
fun main(args: Array<String>) {
    val big = LensAtlas.read(File(args[0]).nameWithoutExtension.removePrefix("atlante_"), File(args[0]).readBytes())
    val fine = LensAtlas.read(File(args[1]).nameWithoutExtension.removePrefix("atlante_"), File(args[1]).readBytes()).also { it.outer = big }
    val lib = LensRadarLibreria.leggi(File(args[2]))
    val maps = args[3]; val kmA = args[5].toDouble(); val kmB = args[6].toDouble()
    val coast = CoastIndex(MapPack.read(File("$maps/europa.rgm").inputStream()), MapPack.read(File("$maps/world.rgm").inputStream()), listOf(MapPack.read(File("$maps/atlantico.rgm").inputStream())))
    WindGrid.landTest = { la, lo -> coast.onLand(la, lo) }          // come l'app: il vento dei GRIB dai soli punti di mare
    val scene = args.drop(7).chunked(4)
    val out = DataOutputStream(BufferedOutputStream(FileOutputStream(args[4])))
    fun wi(x: Int) = out.writeInt(Integer.reverseBytes(x))
    fun wd(x: Double) = out.writeLong(java.lang.Long.reverseBytes(java.lang.Double.doubleToLongBits(x)))
    fun wf(x: Float) = out.writeInt(Integer.reverseBytes(java.lang.Float.floatToIntBits(x)))
    out.write("SCH1".toByteArray()); wi(scene.size); wi(fine.nx); wi(fine.ny); for (b in fine.box) wd(b)
    val nx = fine.nx; val ny = fine.ny
    for ((grib, tS, u925S, v925S) in scene) {
        val t = tS.toLong()
        val g = WindGrid.from(GribReader.read(Unpack.grib(File(grib).readBytes())))
        require(g != null && t in g.startMs..g.endMs) { "il GRIB $grib non copre l'ora" }
        // il primo Lens (aria mescolata, come nell'app con GRIB senza profilo) sull'area grande e sull'area fine
        val hB = big.hour(g, null, t)!!
        val hF = fine.hour(g, null, t)!!
        // Lens 2: la situazione e il fattore della libreria (sulle celle del radar)
        val mu = FloatArray(lib.nc) { Float.NaN }; val mv = FloatArray(lib.nc) { Float.NaN }
        for (c in 0 until lib.nc) {
            val w = g.wind(lib.latCella(c), lib.lonCella(c), t) ?: continue
            val ms = w.twsKn / KNOTS_PER_MS; val a = Math.toRadians(w.twdDeg)
            mu[c] = (-ms * sin(a)).toFloat(); mv[c] = (-ms * cos(a)).toFloat()
        }
        // senza il vento a 925 hPa ("nan"), come Shift: la media del vento del modello nelle celle
        var u925 = u925S.toDoubleOrNull() ?: Double.NaN; var v925 = v925S.toDoubleOrNull() ?: Double.NaN
        if (u925.isNaN() || v925.isNaN()) { val ok = (0 until lib.nc).filter { !mu[it].isNaN() }; u925 = ok.sumOf { mu[it].toDouble() } / ok.size; v925 = ok.sumOf { mv[it].toDouble() } / ok.size }
        val s = LensRadar.situazione(lib, mu, mv, u925, v925, ((t / 3_600_000L) % 24L).toInt())
        val g0 = lib.fattore(s, FloatArray(lib.nc) { LensRadar.provenienza(mu[it], mv[it]) })
        // il modello sull'atlante grande (righe da sud, come Lens2Use.modelOnAtlas)
        val bU = FloatArray(big.nx * big.ny); val bV = FloatArray(big.nx * big.ny)
        for (j in 0 until big.ny) for (i in 0 until big.nx) {
            val q = j * big.nx + i
            val w = g.wind(big.box[0] + (j + 0.5) * big.cellLat, big.box[2] + (i + 0.5) * big.cellLon, t)
            if (w == null) { bU[q] = Float.NaN; bV[q] = Float.NaN; continue }
            val ms = w.twsKn / KNOTS_PER_MS; val a = Math.toRadians(w.twdDeg)
            bU[q] = (-ms * sin(a)).toFloat(); bV[q] = (-ms * cos(a)).toFloat()
        }
        // C: il livello misurato sul Lens dell'area grande, come Shift 0.59.1 (Lens2Use.apply con levelRef)
        val liv30 = LensLivello(lib, g0, big.speedSampler(hB.u, hB.v), big.speedSampler(bU, bV, rowsFromSouth = true), kmA)
        val liv5 = LensLivello(lib, g0, big.speedSampler(hB.u, hB.v), big.speedSampler(bU, bV, rowsFromSouth = true), kmB)
        // B: il radar cella per cella reso relativo al Lens dell'area fine (Shift 0.58)
        val gRel = LensRadar.relativoALens(lib, g0, FloatArray(lib.nc) { hypot(mu[it], mv[it]) }, fine.speedSampler(hF.u, hF.v))
        val campo = LensRadarCampo(lib, gRel)
        val sea = fine.seaInner()
        val fld = Array(11) { FloatArray(nx * ny) }
        val m10 = fine.model10(g, null, t)!!
        for (r in 0 until ny) {
            val lat = fine.box[1] - (r + 0.5) * fine.cellLat
            for (c in 0 until nx) {
                val lon = fine.box[2] + (c + 0.5) * fine.cellLon; val p = r * nx + c
                val w = g.wind(lat, lon, t)
                if (w != null) { val ms = w.twsKn / KNOTS_PER_MS; val a = Math.toRadians(w.twdDeg)
                    fld[0][p] = ms.toFloat(); fld[1][p] = (-ms * sin(a)).toFloat(); fld[2][p] = (-ms * cos(a)).toFloat() } else { fld[0][p] = Float.NaN; fld[1][p] = Float.NaN; fld[2][p] = Float.NaN }
                fld[3][p] = hF.u[p]; fld[4][p] = hF.v[p]
                fld[5][p] = campo.fattore(lat, lon)
                fld[6][p] = liv30.factor(lat, lon).toFloat(); fld[7][p] = liv5.factor(lat, lon).toFloat()
                fld[8][p] = if (sea[p]) 1f else 0f
                fld[9][p] = g0.filter { !it.isNaN() }.average().toFloat()
                fld[10][p] = m10.s[(r + fine.margin) * fine.NX + c + fine.margin].toFloat()
            }
        }
        for (a in fld) for (x in a) wf(x)
        println("scena ${java.time.Instant.ofEpochMilli(t)}: vento al largo ${"%.1f".format(hypot(s.a0, s.b0) * 1.944)} nodi da ${"%.0f".format((Math.toDegrees(atan2(s.a0, s.b0)) + 360) % 360)}°, " +
            "fattore medio della libreria ${"%.2f".format(fld[9][0])}, Fr di A ${"%.2f".format(hF.froude)}")
    }
    out.close()
}
