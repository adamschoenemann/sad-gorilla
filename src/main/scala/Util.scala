
// Some utility classes and definitions


object Util {

  // just a thin wrapper of a Map
  case class LookupTable(map:Map[(Char,Char),Int]) {

    def lookup(i:Char, j:Char):Int = {

      try {
        map.get((i,j)).get
      } catch {
        case e:NoSuchElementException => {
          println(s"Could not find ($i,$j)")
          throw e
        }
      }

    }
  }

  // helper method for printing flat 2d arrays
  private def print2d[A](arr:Seq[A], rows:Int, cols:Int) = {
    for (r <- 0 until rows)
      println(arr.slice(r*cols, r*cols+cols).mkString("\t"))
  }

  // a 2D array, implemented as a flat array with rows and cols
  case class Array2D[A](arr:Array[A], rows:Int, cols:Int) extends Get2D[A] {

    def get(r:Int, c:Int):A =
      try {arr(cols*r + c)}
      catch {
        case (e:IndexOutOfBoundsException) => {
          println("Out of bounds: ", (r,c))
          throw e
        }
      }

    def set(r:Int, c:Int, x:A) =
      arr(r * cols + c) = x

    def print() = print2d(arr, rows, cols)

  }

  // A trait that generalizes 2d things
  trait Get2D[A] {

    def get(r:Int, c:Int):A

  }

  // an implicit conversion from maps to Get2D
  implicit def map2dToGet2D[A](map:Map[(Int,Int),A]):Get2D[A] =
    new Get2D[A] {
      def get(r:Int, c:Int) = map((r,c))
    }

  // a helper method for unfolding stuff to lists! Basically, the opposite
  // of foldRight
  def unfoldRight[A,B] (seed:B) (f: B => Option[(A, B)]):List[A] = {
    f(seed).map({case (a,s2) => a :: unfoldRight(s2)(f)}).getOrElse(Nil)
  }

}
