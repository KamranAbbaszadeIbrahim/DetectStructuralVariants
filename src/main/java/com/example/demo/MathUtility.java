package com.example.demo;

import java.util.List;

/**
 * MathUtility class - class with predefined math statistics functions.
 */
public abstract class MathUtility {
    /**
     * Lanczos Gamma Function approximation - Coefficients
      */
    private static final double[] lgfCoeff = {1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179E-2, -0.5395239384953E-5};

    /**
     * Lanczos Gamma Function approximation - small gamma
     */
    private static final double lgfGamma = 5.0;

    /**
     * Lanczos Gamma Function approximation - N (number of coefficients -1)
     */
    private static final int lgfN = 6;

    /**
     * A small number close to the smallest representable floating point number
     */
    public static final double FPMIN = 1e-300;

    /**
     * tolerance used in the contFract method
     */
    private static final double cfTol = 1.0e-8;

    /**
     * maximum number of iterations allowed in the contFract method
     */
    private static final int cfMaxIter = 500;

    /**
     * meanInt
     *
     * @param list
     * @return Double
     */
    public static Double meanInt(List<Integer> list){
        double mean = 0;
        int count = list.size();
        for(Integer i: list){
            mean += i;
        }
        return mean / count;
    }

    /**
     * meanDouble
     *
     * @param list
     * @return Double
     */
    public static Double meanDouble(List<Double> list){
        double mean = 0;
        int count = list.size();
        for(Double i: list){
            mean += i;
        }
        return mean / count;
    }

    /**
     * standardDeviationInt
     *
     * @param list
     * @return Double
     */
    public static Double standardDeviationInt(List<Integer> list){
        double mean = meanInt(list);
        double standardDeviation = 0;

        for(Integer element: list){
            standardDeviation = standardDeviation + Math.pow(element - mean,2);
        }
        standardDeviation = standardDeviation / list.size();
        return Math.sqrt(standardDeviation);
    }

    /**
     * standardDeviationDouble
     *
     * @param list
     * @return Double
     */
    public static Double standardDeviationDouble(List<Double> list){
        double mean = meanDouble(list);
        double standardDeviation = 0;

        for(Double element: list){
            standardDeviation = standardDeviation + Math.pow(element - mean,2);
        }
        standardDeviation = standardDeviation / list.size();
        return Math.sqrt(standardDeviation);
    }

    /**
     * Returns the Student's t cumulative distribution function probability
     * @param tValue
     * @param df
     * @return double
     */
    public static double studentTCDF(double tValue, double df){
        if(tValue!=tValue)throw new IllegalArgumentException("argument tValue is not a number (NaN)");

        if(tValue==Double.POSITIVE_INFINITY){
            return 1.0;
        }
        else{
            if(tValue==Double.NEGATIVE_INFINITY){
                return 0.0;
            }
            else{
                double ddf = (double)df;
                double x = ddf/(ddf+tValue*tValue);
                return 0.5D*(1.0D + (regularisedBetaFunction(ddf/2.0D, 0.5D, 1) - regularisedBetaFunction(ddf/2.0D, 0.5D, x)) * sign(tValue));
            }
        }
    }

    /**
     * SIGN
     * returns -1 if x < 0 else returns 1
     * double version
     * @param x
     * @return double
     */
    public static double sign(double x){
        if (x<0.0){
            return -1.0;
        }
        else{
            return 1.0;
        }
    }

    /**
     * returns -1 if x < 0 else returns 1
     * int version
     * @param x
     * @return int
     */
    public static int sign(int x){
        if (x<0){
            return -1;
        }
        else{
            return 1;
        }
    }

    /**
     * returns -1 if x < 0 else returns 1
     * long version
     * @param x
     * @return long
     */
    public static long sign(long x){
        if (x<0){
            return -1;
        }
        else{
            return 1;
        }
    }

    /**
     * Regularised Incomplete Beta function
     * Continued Fraction approximation (see Numerical recipies for details of method)
     * @param z
     * @param w
     * @param x
     * @return double
     */
    public static double regularisedBetaFunction(double z, double w, double x){
        if(x<0.0D || x>1.0D)throw new IllegalArgumentException("Argument x, "+x+", must be lie between 0 and 1 (inclusive)");
        double ibeta = 0.0D;
        if(x==0.0D){
            ibeta=0.0D;
        }
        else{
            if(x==1.0D){
                ibeta=1.0D;
            }
            else{
                // Term before continued fraction
                ibeta = Math.exp(logGamma(z+w) - logGamma(z) - logGamma(w) + z*Math.log(x) + w*Math.log(1.0D-x));
                // Continued fraction
                if(x < (z+1.0D)/(z+w+2.0D)){
                    ibeta = ibeta*contFract(z, w, x)/z;
                }
                else{
                    // Use symmetry relationship
                    ibeta = 1.0D - ibeta*contFract(w, z, 1.0D-x)/w;
                }
            }
        }
        return ibeta;
    }

    /**
     * Incomplete fraction summation used in the method regularisedBetaFunction
     * modified Lentz's method
     * @param a
     * @param b
     * @param x
     * @return double
     */
    public static double contFract(double a, double b, double x){

        double aplusb = a + b;
        double aplus1 = a + 1.0D;
        double aminus1 = a - 1.0D;
        double c = 1.0D;
        double d = 1.0D - aplusb*x/aplus1;
        if(Math.abs(d)<FPMIN)d = FPMIN;
        d = 1.0D/d;
        double h = d;
        double aa = 0.0D;
        double del = 0.0D;
        int i=1, i2=0;
        boolean test=true;
        while(test){
            i2=2*i;
            aa = i*(b-i)*x/((aminus1+i2)*(a+i2));
            d = 1.0D + aa*d;
            if(Math.abs(d)<FPMIN)d = FPMIN;
            c = 1.0D + aa/c;
            if(Math.abs(c)<FPMIN)c = FPMIN;
            d = 1.0D/d;
            h *= d*c;
            aa = -(a+i)*(aplusb+i)*x/((a+i2)*(aplus1+i2));
            d = 1.0D + aa*d;
            if(Math.abs(d)<FPMIN)d = FPMIN;
            c = 1.0D + aa/c;
            if(Math.abs(c)<FPMIN)c = FPMIN;
            d = 1.0D/d;
            del = d*c;
            h *= del;
            i++;
            if(Math.abs(del-1.0D) < cfTol)test=false;
            if(i>cfMaxIter){
                test=false;
                System.out.println("Maximum number of iterations ("+cfMaxIter+") exceeded in Stat.contFract in Stat.incompleteBeta");
            }
        }
        return h;

    }


    /**
     * log to base e of the Gamma function
     * Lanczos approximation (6 terms)
     * Retained for backward compatibility
     * @param x
     * @return double
     */
    public static double logGamma(double x){
        double xcopy = x;
        double fg = 0.0D;
        double first = x + lgfGamma + 0.5;
        double second = lgfCoeff[0];

        if(x>=0.0){
            first -= (x + 0.5)*Math.log(first);
            for(int i=1; i<=lgfN; i++)second += lgfCoeff[i]/++xcopy;
            fg = Math.log(Math.sqrt(2.0*Math.PI)*second/x) - first;
        }
        else{
            fg = Math.PI/(gamma(1.0D-x)*Math.sin(Math.PI*x));

            if(fg!=1.0/0.0 && fg!=-1.0/0.0){
                if(fg<0){
                    throw new IllegalArgumentException("\nThe gamma function is negative");
                }
                else{
                    fg = Math.log(fg);
                }
            }
        }
        return fg;
    }

    /**
     * Gamma function
     * Lanczos approximation (6 terms)
     * retained for backward compatibity
     * @param x
     * @return double
     */
    public static double gamma(double x){

        double xcopy = x;
        double first = x + lgfGamma + 0.5;
        double second = lgfCoeff[0];
        double fg = 0.0D;

        if(x>=0.0){
            first = Math.pow(first, x + 0.5)*Math.exp(-first);
            for(int i=1; i<=lgfN; i++)second += lgfCoeff[i]/++xcopy;
            fg = first*Math.sqrt(2.0*Math.PI)*second/x;
        }
        else{
            fg = -Math.PI / (x * gamma(-x) * Math.sin(Math.PI * x));
        }
        return fg;
    }
}
