// TODO: Make it work properly with stdin!

import scala.collection.mutable.HashMap
import io.Source
import java.lang.Math.{max}
import Util._

object Main {

  // just a type alias for clarity
  type Protein = Seq[Char]

  def main(args:Array[String]) {

    val blosumSource = Source.fromFile("./data/BLOSUM62.txt")
    val fastasSrc = getSource(args)
    try {
      val lookup = parseBLOSUM(blosumSource.getLines.toSeq)
      val shiftCost = lookup.lookup('A', '-')
      val (human :: fastas) = parseFASTAs(fastasSrc.getLines)
      val humanProtein = human.proteins
      fastas.map(fasta => {
        fasta.label -> solveFun (lookup, shiftCost) (humanProtein, fasta.proteins)
      }).sortBy(- _._2._1).foreach (result => {
        val (label, (score, al1, al2)) = result
        println(human.label + "--" + label + ": " + score)
        println(al1)
        println(al2)
        println(List.fill (30) ('-') mkString)
      })
    } finally {
      blosumSource.close()
      fastasSrc.close();
    }

  }

  def getSource(args:Array[String]):Source = {
    if (args.length == 1 && args(0) == "test")
      Source.fromFile("./data/HbB_FASTAs-in.txt")
    else
      args.headOption.map(Source.fromFile(_))
          .getOrElse(Source.stdin)
  }

  // a FASTA has a label and a sequence of proteins
  case class FASTA(label:String, proteins:IndexedSeq[Char])

  // parse the fastas from a stream of lines
  def parseFASTAs(iter:Iterator[String]):Seq[FASTA] = {
    iter.grouped(4).foldLeft[List[FASTA]] (Nil) ((acc, fasta) => {
      assert (fasta(0).startsWith(">"))
      val label = fasta.head.tail.takeWhile(!_.isDigit).trim
      val proteins = fasta.tail.mkString.toIndexedSeq
      FASTA(label, proteins) :: acc
    }).reverse
  }


  // imperative implementation of assemble, which assembles a solution from
  // the "memory" generated from a solving function
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
  def assembleFun(table:LookupTable, memory:Get2D[Int], delta:Int)
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

  // an imperative bottom-up implementation of the string-matching algorithm
  def solve(table:LookupTable, delta:Int)(p1: Protein, p2:Protein):(Int,String,String) = {

    val rows = p1.length + 1
    val cols = p2.length + 1
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

  }

  // A recursive, functional implementation of the string-matching algorithm
  def solveFun(table:LookupTable, delta:Int)(p1:Protein, p2:Protein): (Int,String,String) = {
    type Memory = Map[(Int,Int),Int]

    // import home-rolled State monad :D
    import State._
    def helper(xs:Protein, ys:Protein):State[Memory, Int] = {

      def help(memory:Memory):State[Memory,Int] = {
        if (xs.isEmpty) {
          val r = ys.length * delta
          for {
            _ <- put (memory.updated((0, ys.length), r))
          } yield r
        } else if (ys.isEmpty) {
          val r = xs.length * delta
          for {
            _ <- put (memory.updated((xs.length, 0), r))
          } yield r
        } else if (memory.contains((xs.length, ys.length))) {
          unit(memory.get((xs.length, ys.length)).get)
        } else {
          val x = xs.last
          val y = ys.last
          val alpha = table.lookup(x, y)
          for {
            r1 <- helper (xs.init, ys.init)
            r2 <- helper (xs.init, ys)
            r3 <- helper (xs, ys.init)
            memory2 <- get[Memory] // IMPORTANT. get an "updated" version of memory here
            m = max(alpha + r1, max(delta + r2, delta + r3))
            _ <- put (memory2.updated((xs.length, ys.length), m))
          } yield m
        }
      }

      for {
        memory <- get
        m <- help (memory)
      } yield m

    }

    val (score, memory) = helper(p1,p2).run(Map[(Int,Int),Int]())
    val (s1, s2) = assembleFun(table, memory, delta) (p1, p2)
    (score, s1, s2)

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