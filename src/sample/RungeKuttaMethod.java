package sample;
import javafx.scene.chart.XYChart;

public class RungeKuttaMethod extends NumericalMethods {

    private double RungeKuttaEq(double x, double h, double y) {
        double k1 = dydx(x,y);
        double k2 = dydx((x + h/2),(y + k1/2));
        double k3 = dydx((x + h/2),(y + k2/2));
        double k4 = dydx((x + h),(y + k3));

        if (k1 == CONST || k2 == CONST || k3 == CONST || k4 == CONST) { return CONST; }
        else { return (y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)); }
    }

    @Override
    XYChart.Series[] getGraphic(double x0, double y0, double x1, int N) {
        double tmp = x0;
        double h = ((x1 - x0) / N);

        double[] yEx;
        double[] yRK = new double[N + 1];
        double[] globErr = new double[N + 1];
        double[] localErr = new double[N + 1];
        yRK[0] = y0;

        XYChart.Series ESn = new XYChart.Series();  ESn.setName("Exact");         // Exact solution with -x
        XYChart.Series ESp = new XYChart.Series();  ESp.setName("Exact");         // Exact solution with +x
        XYChart.Series RKM = new XYChart.Series();   RKM.setName("Runge-Kutta");  // Euler's numerical method
        XYChart.Series GE = new XYChart.Series();   GE.setName("Global");         // Global Error with
        XYChart.Series LE = new XYChart.Series();   LE.setName("Local");          // Local Error with positive and negative are same

        RKM.getData().add(new XYChart.Data<>(x0, yRK[0]));

        x0+=h;
        int trap = -1;
        for(int i=1; i<=N; i++, x0+=h){
            x0 = Math.round(x0 * (1/h)) / (1/h);
            yRK[i] = RungeKuttaEq(x0, h, yRK[i-1]);
            if(yRK[i] == CONST) {trap = i; break;}
            if(yRK[i] > BOUND*2 || yRK[i] < -BOUND*2) continue;
            RKM.getData().add(new XYChart.Data<>(x0, yRK[i]));
        }
        XYChart.Series[] res = {RKM, RKM, RKM, RKM, RKM};

        if(y0*y0 > 2) {
            x0 = tmp;
            yEx = ExactSolution(x0,y0,x1,N);

            globErr[0] = Math.abs((yEx[0] - yRK[0]));
            GE.getData().add(new XYChart.Data<>(x0, globErr[0]));
            x0+=h;
            for(int i=1; i<=N; i++, x0+=h) {
                if(trap == i) break;
                if(yEx[i] > BOUND || yRK[i] > BOUND || yRK[i]< -BOUND) continue;

                globErr[i] = Math.abs((yEx[i] - yRK[i]));
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
                yNS[j] = RungeKuttaEq(x0, h, yNS[j-1]);
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
