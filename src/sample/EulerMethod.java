package sample;
import javafx.scene.chart.XYChart;

public class EulerMethod extends NumericalMethods {

    private double EulerEq(double x, double h, double y) {
        if (x== 0 || y == 0) { return CONST; }
        else { return (y + h*dydx(x,y)); }
    }

    @Override
    XYChart.Series[] getGraphic(double x0, double y0, double x1, int N) {
        double tmp = x0;
        double h = ((x1 - x0) / N);

        double[] yEx;
        double[] yEul = new double[N + 1];
        double[] globErr = new double[N + 1];
        double[] localErr = new double[N + 1];
        yEul[0] = y0;

        XYChart.Series ESn = new XYChart.Series();  ESn.setName("Exact");    // positive Exact solution
        XYChart.Series ESp = new XYChart.Series();  ESp.setName("Exact");    // negative Exact solution
        XYChart.Series EM = new XYChart.Series();   EM.setName("Euler's");   // Euler's numerical method
        XYChart.Series GE = new XYChart.Series();   GE.setName("Global");    // Global Error with negative Exact solution
        XYChart.Series LE = new XYChart.Series();   LE.setName("Local");     // Local Error with positive and negative are same

        EM.getData().add(new XYChart.Data<>(x0, yEul[0]));

        x0+=h;
        int trap = -1;
        for(int i=1; i<=N; i++, x0+=h){
            x0 = Math.round(x0 * (1/h)) / (1/h);
            yEul[i] = EulerEq(x0, h, yEul[i-1]);
            if(yEul[i] == CONST) {trap = i; break;}
            EM.getData().add(new XYChart.Data<>(x0, yEul[i]));
        }

        XYChart.Series[] res = {EM, EM, EM, EM, EM};

        if(y0*y0 > 2) {
            x0 = tmp;
            yEx = ExactSolution(x0,y0,x1,N);

            globErr[0] = Math.abs((yEx[0] - yEul[0]));
            GE.getData().add(new XYChart.Data<>(x0, globErr[0]));

            x0+=h;
            for(int i=1; i<=N; i++, x0+=h) {
                if(trap == i) break;
                if(yEx[i] > BOUND || yEul[i] > BOUND) continue;

                globErr[i] = Math.abs((yEx[i] - yEul[i]));
                localErr[i] = Math.abs(globErr[i] - globErr[i - 1]);
                GE.getData().add(new XYChart.Data<>(x0, globErr[i]));
                LE.getData().add(new XYChart.Data<>(x0, localErr[i]));
            }

            x0 = tmp;
            boolean trapped = false;
            for(int i=0; i<=N; i++, x0+=h) {
                if(yEx[i] == CONST) {trapped = true; continue;}
                if(yEx[i] > BOUND) continue;

                if(!trapped) ESn.getData().add(new XYChart.Data<>(x0, yEx[i]));
                else ESp.getData().add(new XYChart.Data<>(x0, yEx[i]));
            }

            res[0] = ESn;   res[1] = ESn;   res[3] = GE;    res[4] = LE;    if (trapped) res[1] = ESp;
        }
        return res;
    }

    @Override
    XYChart.Series TotalError(double x0, double y0, double x1, int Nmin, int Nmax) {
        XYChart.Series TE = new XYChart.Series();  TE.setName("Total Error");
        double tmp = x0;

        for(int i = Nmin, no = 0; i<=Nmax; i++, no++) {
            x0 = tmp;
            double h = ((x1 - x0) / i);
            double[] yEx =  ExactSolution(x0, y0, x1, i);;
            double[] yNS = new double[i+1];
            double[] GE = new double [i+1];
            yNS[0] = y0;

            int trap = -1;  x0+=h;
            for(int j=1; j<=i; j++, x0+=h) {
                x0 = Math.round(x0 * (1/h)) / (1/h);
                yNS[j] = EulerEq(x0, h, yNS[j-1]);
                if(yNS[i] == CONST) {trap = j; break;}
            }

            GE[0] = Math.abs(yEx[0] - yNS[0]);
            double max = -999999, check;
            x0 = tmp+h;
            for(int j=1; j<=i; j++, x0+=h){
                if(trap == j) break;
                GE[j] = Math.abs(yEx[j] - yNS[j]);
                check =  Math.abs(GE[j] - GE[j-1]);
                if (check > max) max = check;
            }
            TE.getData().add(new XYChart.Data<>(no, max));
        }
        return TE;
    }
}
