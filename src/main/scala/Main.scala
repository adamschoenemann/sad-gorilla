
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import io.Source
import util.matching._
import java.lang.Math.{max}


object Main {

  type Char2dMap = Map[(Char,Char),Int]

  case class LookupTable(map:Char2dMap) {

    def lookup(i:Char, j:Char):Int = {

      map.get((i,j)).get
    }

  }

  type Protein = Seq[Char]

  case class Array2D[A](arr:Array[A], rows:Int, cols:Int) {

    // private val arr:Array[A] = Array.fill(cols*rows)(default)

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

    private def print2d(arr:Array[A], rows:Int, cols:Int) = {
      for (r <- 0 until rows)
        println(arr.slice(r*cols, r*cols+cols).mkString("\t"))
    }
  }

  def testLookup = {
    val map = new HashMap[(Char,Char), Int]
    map.put(('A','A'),2)
    map.put(('A','C'),-1)
    map.put(('A','G'),1)
    map.put(('A','T'),-1)
    map.put(('C','A'),-1)
    map.put(('C','C'),2)
    map.put(('C','G'),-1)
    map.put(('C','T'),1)
    map.put(('G','A'),1)
    map.put(('G','C'),-1)
    map.put(('G','G'),2)
    map.put(('G','T'),-1)
    map.put(('T','A'),-1)
    map.put(('T','C'),1)
    map.put(('T','G'),-1)
    map.put(('T','T'),2)
    LookupTable(map.toMap)
  }

  def main(args:Array[String]) {

    val source = Source.fromFile("./data/BLOSUM62.txt")
    try {
      val lookup = parseBLOSUM(source.getLines.toSeq)
      val sphinx = "KQRK".toList
      val snark = "KQRIKAAKABK".toList
      val s = "ACGGTAG".toIndexedSeq
      val t = "CCTAAG".toIndexedSeq
      val r = solve (lookup, -4) (sphinx, snark)
      println(r)
    } finally source.close()

  }

  def unfoldRight[A,B] (seed:B) (f: B => Option[(A, B)]):List[A] = {
    f(seed).map({case (a,s2) => a :: unfoldRight(s2)(f)}).getOrElse(Nil)
  }

  def assembleImp(table:LookupTable, memory:Array2D[Int], delta:Int)
                 (p1:Protein, p2:Protein)
                 :(String,String) = {
    var i = p1.length
    var j = p2.length
    var s1 = List[Char]()
    var s2 = List[Char]()
    while (i > 0 && j > 0) {
      if (memory.get(i,j) == table.lookup(p1(i-1), p2(j-1)) + memory.get(i-1,j-1)) {
        s1 = p1(i-1) :: s1
        s2 = p2(j-1) :: s2
        i = i - 1
        j = j - 1
      } else if (memory.get(i,j) == delta + memory.get(i-1,j)) {
        s1 = p1(i-1) :: s1
        s2 = '-' :: s2
        i = i - 1
      } else {
        s1 = '-' :: s1
        s2 = p2(j-1) :: s2
        j = j - 1
      }

    }
    (s1.mkString,s2.mkString)
  }

  // functional implementation of assemble
  def assembleFun(table:LookupTable, memory:Array2D[Int], delta:Int)
                 (p1:Protein, p2:Protein)
                 :(String,String) = {

    val rows = p1.length
    val cols = p2.length
    val (s1,s2) = unfoldRight[(Char,Char), (Int,Int)]((rows, cols)) {seed => {
      val (i, j) = seed
      if (i > 0 && j > 0) {
        if (memory.get(i,j) == table.lookup(p1(i-1), p2(j-1)) + memory.get(i-1,j-1)) {
          val a = (p1(i-1), p2(j-1))
          val seed2 = (i-1, j-1)
          Some((a,seed2))
        } else if (memory.get(i,j) == delta + memory.get(i-1,j)) {
          val a = (p1(i-1), '-')
          val seed2 = (i-1, j)
          Some((a,seed2))
        } else {
          val a = ('-', p2(j-1))
          val seed2 = (i, j-1)
          Some((a,seed2))
        }
      } else None
    }}.foldLeft[(List[Char],List[Char])] ((Nil,Nil)) ((acc, x) => x match {
      case (a,b) => (a :: acc._1, b :: acc._2)
    })

    (s1.mkString, s2.mkString)
  }

  def solve(table:LookupTable, delta:Int)(p1: Protein, p2:Protein):(Int,String,String) = {

    val rows = p1.length + 1
    val cols = p2.length + 1
    // val memory = Array.fill (cols*rows) (0)
    val memory = Array2D(Array.fill(rows*cols)(0), rows, cols)


    for (i <- 0 until rows)
      memory.set(i, 0, delta * i)
    for (i <- 0 until cols)
      memory.set(0, i, delta * i)

    for (i <- 1 until rows) {
      for (j <- 1 until cols) {
        val ismatch = table.lookup(p1(i-1), p2(j-1)) + memory.get(i-1,j-1)

        val dash1 = delta + memory.get(i  , j-1)
        val dash2 = delta + memory.get(i-1, j)
        val m     = max(ismatch, max(dash1, dash2))
        memory.set(i,j, m)
      }
    }

    val score = memory.get(rows-1, cols-1)

    val xs = unfoldRight(0) (x => if (x < 10) Some((x,x+1)) else None)

    val (s1, s2) = assembleImp(table, memory, delta) (p1, p2)

    (score, s1, s2)

    // ------------------------------------------------
    // Recursive impl!
    // ------------------------------------------------
    // def helper(xs:Protein, ys:Protein):Int = {
    //   if (memory.contains((xs.length, ys.length)))
    //     memory.get((xs.length, ys.length)).get
    //   else {
    //     val result =
    //       if (xs.isEmpty) {
    //         // inverse.put((xs.length, ys.length), ('-',ys.lastOption.getOrElse('-')))
    //         ys.length * (-4)
    //       } else if (ys.isEmpty) {
    //         // inverse.put((xs.length, ys.length), (xs.lastOption.getOrElse('-'),'-'))
    //         xs.length * (-4)
    //       } else {
    //         val x = xs.last
    //         val y = ys.last
    //         val alpha = table.lookup(x, y)
    //         val del = -4
    //         val r = List(alpha + helper(xs.init, ys.init),
    //           del + helper(xs.init, ys),
    //           del + helper(xs, ys.init)
    //         ).max
    //         // inverse.put((xs.length, ys.length), (x, y))
    //         r
    //       }
    //     memory.put((xs.length, ys.length), result)

    //     result
    //   }

    // }

    // var out1 = ""
    // var out2 = ""


    // def findsol(i:Int, j:Int) = {
    //   val a = memory((i,j))
    //   val alpha = table.lookup(p1(i), p2(j))
    //   if (a == alpha + memory((i-1, j-1))) {
    //     out1 = out1 + p1(i)
    //     out2 = out2 + p2(j)
    //   }
    // }


    // val r = helper(p1, p2)
    // // println(inverse((p1.length-2, p2.length-2)))
    // // println(inverse)
    // for (i <- 1 until p1.length) {
    //   for (j <- 1 until p2.length) {
    //     findsol(i,j)
    //   }
    // }
    // println(out1)
    // r
  }

  def mkSolver(table:LookupTable) (memory:HashMap[(Char,Char),Int])
           :(Protein, Protein) => (Int, (String,String)) = {

    // def helper(xs:Protein, ys:Protein):(Int, (String,String)) = {
    //   if (xs.isEmpty)
    //     (ys.length * (-4), (List.fill(ys.length)('-').mkString, ys.mkString))
    //   else if (ys.isEmpty)
    //     (xs.length * (-4), (xs.mkString, List.fill(xs.length)('-').mkString))
    //   else {
    //     val x = xs.last
    //     val y = ys.last
    //     val combs = List(('-',y),(x,y),(x,'-'))
    //       .map({case (a,b) => {
    //         val (m1, v1) = helper (xs.init, ys.init)
    //         (table.lookup(a,b), (a, b))
    //       }})
    //     val max@(m,(a,b)) = combs.maxBy(_._1)
    //     memory.put((a,b),m)
    //     (m, (a.toString, b.toString))
    //   }
    // }

    // helper
    // val m = Math.min(table.lookup('-',y), table.lookup(x,y), table.lookup(x,'-'))
    ???
  }

  def parseBLOSUM(lines:Seq[String]):LookupTable = {
    def parseLine(line:String) =
      "\\s+".r.replaceAllIn(line.trim, ":")
              .replace("*", "-").split(':').toIndexedSeq

    val filtered = lines.filter(! _.trim.startsWith("#"))

    val first = filtered.head
    val rest  = filtered.tail

    val header = parseLine(first).map(_.charAt(0))
    // println(header)
    val map = new HashMap[(Char,Char),Int]()
    for (line <- rest) {
      val vals = parseLine(line)
      val row = vals.head.charAt(0)
      for ((v,i) <- vals.tail.zipWithIndex) {
        val col = header(i)
        map.put((row,col), v.toInt)
      }
      // println(vals)
    }

    // println(map)
    LookupTable(map.toMap)

  }


}