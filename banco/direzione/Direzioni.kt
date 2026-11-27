import app.regata.core.*
import java.io.*
import kotlin.math.*

/**
 * La direzione di Lens nei punti delle stazioni (27/11 sera, Giano: "sono le rotazioni quello che contano di più").
 * Uso: Direzioni <atlante> <cartella mappe> <uscita .txt> <punti: nome:lat:lon,...> poi per scena <GRIB vento> <GRIB profilo> <tempo ms>
 * Per ogni scena e punto scrive una riga: tempo nome Fr_aria Fr_mescolata e, nella cella della stazione e nella cella di
 * mare aperto più vicina (mare con tutto mare attorno per 3 celle), u e v (m/s, verso cui soffia) di:
 *   M   il modello come lo riceve Lens (lisciato su 4 km)
 *   L0  Lens Fast senza la rotazione del rilievo (TERRAIN_DIR_WEIGHT 0)
 *   LA  Lens Fast con la rotazione piena e l'aria della scena (profilo), come il Detailed sceglie la stabilità
 *   LM  Lens Fast con la rotazione piena e l'aria mescolata (GRIB senza profilo)
 *   LS  Lens Fast con la rotazione piena e l'aria molto stabile (N = 0,02)
 */
fun main(args: Array<String>) {
    val atl = LensAtlas.read(File(args[0]).nameWithoutExtension.removePrefix("atlante_"), File(args[0]).readBytes())
    val maps = args[1]
    val coast = CoastIndex(MapPack.read(File("$maps/europa.rgm").inputStream()), MapPack.read(File("$maps/world.rgm").inputStream()), listOf(MapPack.read(File("$maps/atlantico.rgm").inputStream())))
    WindGrid.landTest = { la, lo -> coast.onLand(la, lo) }
    val nx = atl.nx; val ny = atl.ny; val sea = atl.seaInner()
    val aperto = BooleanArray(nx * ny) { p ->
        val r = p / nx; val c = p % nx
        var ok = r >= 3 && r < ny - 3 && c >= 3 && c < nx - 3
        if (ok) for (dr in -3..3) for (dc in -3..3) if (!sea[(r + dr) * nx + c + dc]) ok = false
        ok
    }
    class Punto(val nome: String, val cella: Int, val mare: Int, val kmMare: Double)
    val punti = args[3].split(",").mapNotNull { s ->
        val (nome, la, lo) = s.split(":"); val lat = la.toDouble(); val lon = lo.toDouble()
        val r = ((atl.box[1] - lat) / atl.cellLat).toInt(); val c = ((lon - atl.box[2]) / atl.cellLon).toInt()
        if (r < 0 || r >= ny || c < 0 || c >= nx) return@mapNotNull null
        val kx = cos(Math.toRadians(lat)); var best = -1; var bd = Double.MAX_VALUE
        for (p in 0 until nx * ny) if (aperto[p]) {
            val d = hypot((p / nx - r) * atl.cellLat * 111.19, (p % nx - c) * atl.cellLon * 111.19 * kx)
            if (d < bd) { bd = d; best = p }
        }
        Punto(nome, r * nx + c, best, bd)
    }
    val out = PrintWriter(File(args[2]))
    out.println("# punti: " + punti.joinToString(" ") { "${it.nome}:${"%.2f".format(it.kmMare)}km" })
    var ultimo = ""; var g: WindGrid? = null
    for (sc in args.drop(4).chunked(3)) {
        val (gW, gP, tS) = sc; val t = tS.toLong()
        if (gW != ultimo) { g = WindGrid.from(GribReader.read(Unpack.grib(File(gW).readBytes()))); ultimo = gW }
        val gg = g ?: continue
        if (t !in gg.startMs..gg.endMs) continue
        val air = if (File(gP).exists()) aria(File(gP), gg, atl, t) else null
        val m = atl.model10(gg, null, t) ?: continue
        val vecchio = LensAtlas.TERRAIN_DIR_WEIGHT
        LensAtlas.TERRAIN_DIR_WEIGHT = 0.0
        val h0 = atl.hour(gg, null, t, air?.n ?: LensAtlas.N_MIXED, air)!!
        LensAtlas.TERRAIN_DIR_WEIGHT = 1.0
        val hA = atl.hour(gg, null, t, air?.n ?: LensAtlas.N_MIXED, air)!!
        val hM = atl.hour(gg, null, t)!!
        val hS = atl.hour(gg, null, t, 0.02)!!
        LensAtlas.TERRAIN_DIR_WEIGHT = vecchio
        for (p in punti) {
            val sb = StringBuilder("$t ${p.nome} ${"%.3f".format(air?.fr ?: Double.NaN)} ${"%.3f".format(hM.froude)}")
            for (q in intArrayOf(p.cella, p.mare)) {
                if (q < 0) { repeat(10) { sb.append(" nan") }; continue }
                val qq = (q / nx + atl.margin) * atl.NX + q % nx + atl.margin
                sb.append(" %.3f %.3f".format(m.u[qq], m.v[qq]))
                for (h in listOf(h0, hA, hM, hS)) sb.append(" %.3f %.3f".format(h.u[q], h.v[q]))
            }
            out.println(sb)
        }
    }
    out.close()
}
