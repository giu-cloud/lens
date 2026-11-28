import app.regata.core.*
import java.io.*

/**
 * La maschera di terra di un riquadro dalla costa di Shift (per le strisce del radar).
 * Uso: MascheraTerra <cartella mappe> <S> <N> <O> <E> <passo in gradi> <uscita>
 * Scrive ny, nx (int32) e poi ny*nx byte (1 = terra), righe da sud, colonne da ovest; il punto (j, i) è S + (j+0,5) passo, O + (i+0,5) passo.
 */
fun main(a: Array<String>) {
    val maps = a[0]
    val coast = CoastIndex(MapPack.read(File("$maps/europa.rgm").inputStream()), MapPack.read(File("$maps/world.rgm").inputStream()), listOf(MapPack.read(File("$maps/atlantico.rgm").inputStream())))
    val s = a[1].toDouble(); val n = a[2].toDouble(); val o = a[3].toDouble(); val e = a[4].toDouble(); val p = a[5].toDouble()
    val ny = ((n - s) / p).toInt(); val nx = ((e - o) / p).toInt()
    DataOutputStream(BufferedOutputStream(FileOutputStream(a[6]))).use { out ->
        out.writeInt(ny); out.writeInt(nx)
        for (j in 0 until ny) for (i in 0 until nx) out.writeByte(if (coast.onLand(s + (j + 0.5) * p, o + (i + 0.5) * p)) 1 else 0)
    }
}
