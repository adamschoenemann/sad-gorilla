
object State {

  def unit[S,A](a:A):State[S,A] =
    State(s => (a,s))

  def get[S]:State[S,S] =
    State(s => (s,s))

  def put[S](a:S):State[S,Unit] =
    State(s => ((), a))

  def update[S](fn: S => S):State[S,Unit] =
    State(s => ((), fn(s)))

  // I am not actually using sequence for anything, it was just fun to write
  // two implementations of it.
  // This is the hand-rolled version with explicit recursion
  def sequence[S,A](comps:Seq[State[S,A]]):State[S, Seq[A]] = {

    def helper(s:S, comps:Seq[State[S,A]]):(List[A],S) = comps match {
      case Nil => (Nil, s)
      case (x :: xs) => {
        val (y, s2) = x.run(s)
        val (ys, s3) = helper(s2, xs)
        (y :: ys, s3)
      }
    }

    State(s => helper(s, comps))

  }

  // This is the version that simply uses foldLeft and reverse
  def sequence2[S,A](comps:Seq[State[S,A]]):State[S, Seq[A]] =
    State(s => {
      val (ss, r) = comps.foldLeft[(List[A], S)] ((Nil, s)) ((acc, x) => {
        val (lst, s1) = acc
        val (y, s2) = x.run (s1)
        (y :: lst,s2)
      })
      (ss.reverse, r)
    })

}

case class State[S, A](run: S => (A,S)) {

  def flatMap[B](f: A => State[S,B]):State[S,B] = {
    State(s1 => {
      val (y, s2) = this.run(s1)
      f(y).run(s2)
    })
  }

  def map[B](f: A => B):State[S,B] = {
    State(s1 => {
      val (x, s2) = this.run(s1)
      (f(x), s2)
    })
  }

}