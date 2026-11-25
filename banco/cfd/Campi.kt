import app.regata.core.*
import java.io.File
import java.io.DataOutputStream
import java.io.BufferedOutputStream

// Campi di prova su un atlante: Lens 0.56 di oggi e il primo Lens (fisica piena), per un vento uniforme.
fun main(args: Array<String>) {
    val atl = LensAtlas.read("prova", File(args[0]).readBytes())
    val out = DataOutputStream(BufferedOutputStream(File(args[1]).outputStream()))
    out.writeInt(atl.nx); out.writeInt(atl.ny)
    for (v in atl.box) out.writeDouble(v)
    var k = 2
    while (k + 1 < args.size) {
        val da = args[k].toDouble(); val kn = args[k + 1].toDouble(); k += 2
        val w = object : WindProvider { override val startMs = 0L; override val endMs = Long.MAX_VALUE; override fun wind(lat: Double, lon: Double, timeMs: Long) = Wind(da, kn) }
        val fr = System.getProperty("fr")?.toDoubleOrNull()
        val air = if (fr == null) null else LensAirHour(fr, 0.01, 0.01, kn / 1.944, Double.NaN, Double.NaN)
        val oggi = LensAtlas.physics054 { atl.hour(w, null, 0L, LensAtlas.N_MIXED, air)!! }   // 0.56
        val primo = atl.hour(w, null, 0L, LensAtlas.N_MIXED, air)!!                          // 0.57: il default
        for (h in listOf(oggi, primo)) { for (x in h.u) out.writeFloat(x); for (x in h.v) out.writeFloat(x) }
        println("vento da $da a $kn nodi: Fr ${"%.2f".format(oggi.froude)}")
    }
    out.close()
}
