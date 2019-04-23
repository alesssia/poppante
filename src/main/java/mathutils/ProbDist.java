
package mathutils;

/**
 *
 * ProbDist
 *
 * @author  HkS
 * @version 1.0 beta 2005.07.06
 * @since   1.0
 *
 */

public abstract class ProbDist
{

    public static final int      ITMAX        =     300;
    public static final double   FPMIN        = 1.0E-50;
    public static final double   EPSILON      = 1.0E-20;
//    public static final double   DELTA        = 1.0E-13;
    public static final double   LANCZ_CUTOFF =   700.0;
//    public static final int      MAXPOINTS    =    8000;
    private static final long serialVersionUID = 1L;

  // Costruttori ---------------------------------------------------------------

  /**
   *
   * La classe e' <tt>abstract</tt>, non puo' essere istanziata.
   *
   * @since 1.0
   *
   */

  private ProbDist()
  {
  }

  // Metodo gaussian -----------------------------------------------------------

  /**
   *
   * gaussian
   *
   * @since  1.0
   *
   */

  /* --------------------------------------------------------------------
           double gammq(double a, double x)
  -------------------------------------------------------------------- */
  public static double gammq(double a, double x) {
       double retval, g;

       if (x < 0 || a <= 0)
           retval = -1.0;

       else {
           if (a >= LANCZ_CUTOFF)
               g = LogStirling(a);
           else
               g = Lanczos(a);

           if (x >= a + 1.0)
               retval = gcf(a, x, g);
           else if ((g = gser(a, x, g)) < 0)
               retval = g;
           else
               retval = 1.0 - g;
       }
       return (retval);
  }


  /* --------------------------------------------------------------------
     LogStirling  =  LOG( sqr(2Pi)n^(n+.5)exp(-n) )
  -------------------------------------------------------------------- */
  protected static double LogStirling(double p) {
      return ( 0.5 * Math.log(2.0 * Math.PI) + (0.5 + p) * Math.log(p) - p );
  }

  /* --------------------------------------------------------------
      Lanczos(p)  Lanczos's neat approximation of log(Gamma(p))
                  --slightly modified for speed
  -------------------------------------------------------------- */
  protected static double Lanczos(double p) {
     int j;
     double x, tmp, ser;

     x = p;
     tmp = x + 5.5;
     tmp = tmp - (x + .5) * Math.log(tmp);
     ser = 1.000000000190015 + 76.18009172947146 / ( p + 1.0 );
     ser -= 86.50532032941678 / ( p + 2.0 );
     ser += 24.01409824083091 / ( p + 3.0 );
     ser -= 1.231739572450155 / ( p + 4.0 );
     ser += .001208650973866179 / ( p + 5.0 );
     ser -= 5.395239384953E-06 / ( p + 6.0 );
     return (Math.log(2.506628274631001 * ser / x) - tmp);
  }
  /* --------------------------------------------------------------------
           double gser (a, x, g)
          '
          '  returns the incomplete gamma function P(a,x)
          '  g is Lanczos(a)
          '
  -------------------------------------------------------------------- */
  protected static double gser(double a, double x, double g) {
      int i;
      double sum, del, ap, retval=-1.0;

      if (x == 0.0) {
           retval = 0.0;
      }
      else if (x > 0) {
          ap = a;
          sum = 1.0 / a;
          del = sum;
          for (i=1; i<=ITMAX; i++) {
              ap += 1.0;
              del *= x / ap;
              sum += del;
              if (Math.abs(del) < Math.abs(sum) * EPSILON) {
                  retval = sum * Math.exp(a * Math.log(x) - x - g);
                  break;
              }
         }
      }
      return (retval);
  }
  /* --------------------------------------------------------------------
           double gcf (a, x, g)
     '
     '  returns the incomplete gamma function Q(a,x)
     '  g is Lanczos(a)
     '
  -------------------------------------------------------------------- */
  protected static double gcf(double a, double x, double g) {
      int i;
      double retval, an, b, c, d, del, h;

      b = x + 1.0 - a;
      c = 1.0 / FPMIN;
      d = 1.0 / b;
      h = d;
      for (i=1; i<=ITMAX; i++) {
          an = i * (a - i);
          b += 2.0;
          d = an * d + b;
          if (Math.abs(d) < FPMIN) d = FPMIN;
          c = b + an / c;
          if (Math.abs(c) < FPMIN) c = FPMIN;
          d = 1.0 / d;
          del = d * c;
          h *= del;
          if (Math.abs(del - 1.0) < EPSILON) break;
      }

      if (i > ITMAX) retval = -1.0;
      else retval = Math.exp(a * Math.log(x) - x - g) * h;

      return (retval);
  }


}
