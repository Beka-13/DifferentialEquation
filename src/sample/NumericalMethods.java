package sample;
import javafx.scene.chart.XYChart;

public abstract class NumericalMethods {
    public static final double CONST = -90476.11375;
    public static final double BOUND = 10;

    double dydx(double x, double y) {
        if ( x == 0 || y == 0) { return CONST; }
        else { return ((2 - y*y) / (2*x*x*y)); }
    }

    double[] ExactSolution(double x0, double y0, double x1, int N) {
        double[] yEx = new double[N + 1];
        double c = Math.log(y0*y0 - 2) - (1 / x0);
        double h = ((x1 - x0) / N);

        for (int i=0; i<=N; i++, x0+=h) {
            x0 = Math.round(x0 * (1/h)) / (1/h);
            if(x0 == 0.0) { yEx[i] = CONST; continue; }
            yEx[i] =  Math.sqrt(2 + Math.exp(c + 1/x0));
        }
        return yEx;
    }

    abstract  XYChart.Series[] getGraphic(double x0, double y0, double x1, int N);
    abstract  XYChart.Series TotalError(double x0, double y0, double x1, int Nmin, int Nmax);

}
