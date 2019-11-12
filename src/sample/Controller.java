package sample;

import java.net.URL;
import java.util.ResourceBundle;
import java.util.concurrent.atomic.AtomicReference;

import javafx.fxml.FXML;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.scene.control.Button;
import javafx.event.ActionEvent;

public class Controller {
    @FXML private ResourceBundle resources;
    @FXML private URL location;
    @FXML private LineChart<?, ?> MyChart;
    @FXML private NumberAxis Xaxis;
    @FXML private NumberAxis Yaxis;
    @FXML private LineChart<?, ?> ErrorChart;
    @FXML private NumberAxis ErrorX;
    @FXML private NumberAxis ErrorY;
    @FXML private TextField YoSetter;
    @FXML private TextField XoSetter;
    @FXML private TextField XSetter;
    @FXML private TextField NSetter;
    @FXML private Button SetBtn;
    @FXML private CheckBox EulerChk;
    @FXML private CheckBox ImEulChk;
    @FXML private CheckBox RunKutChk;
    @FXML private LineChart<?, ?> TotalErrorChart;
    @FXML private NumberAxis ErrorX1;
    @FXML private NumberAxis ErrorY1;
    @FXML private TextField NminSetter;
    @FXML private TextField NmaxSetter;
    @FXML private CheckBox Nchk;



    void addSeries(XYChart.Series[] res){
        MyChart.getData().addAll(res[0], res[2]);
        if ( res[0] != res[1]) MyChart.getData().addAll(res[1]);
        ErrorChart.getData().addAll(res[3], res[4]);
    }

    void getEulerData(double x0, double y0, double x1, int N){
        NumericalMethods tmp = new EulerMethod();
        XYChart.Series[] res = tmp.getGraphic(x0,y0,x1,N);
        if(y0*y0 > 2) addSeries(res);
        else MyChart.getData().addAll(res[2]);

        if ( Nchk.isSelected() && (y0*y0 > 2)) {
            int min = Integer.parseInt(NminSetter.getText());
            int max = Integer.parseInt(NmaxSetter.getText());

            XYChart.Series err = tmp.TotalError(x0, y0, x1, min, max);
            TotalErrorChart.getData().addAll(err);
        }
    }

    void getImprovedEulerData(double x0, double y0, double x1, int N){
        NumericalMethods tmp = new ImprovedEuler();
        XYChart.Series[] res = tmp.getGraphic(x0,y0,x1,N);
        if(y0*y0 > 2) addSeries(res);
        else MyChart.getData().addAll(res[2]);

        if ( Nchk.isSelected() && (y0*y0 > 2)) {
            int min = Integer.parseInt(NminSetter.getText());
            int max = Integer.parseInt(NmaxSetter.getText());

            XYChart.Series err = tmp.TotalError(x0, y0, x1, min, max);
            TotalErrorChart.getData().addAll(err);
        }
    }

    void getRungeKuttaData(double x0, double y0, double x1, int N){
        NumericalMethods tmp = new RungeKuttaMethod();
        XYChart.Series[] res = tmp.getGraphic(x0,y0,x1,N);
        if(y0*y0 > 2) addSeries(res);
        else MyChart.getData().addAll(res[2]);

        if ( Nchk.isSelected() && (y0*y0 > 2)) {
            int min = Integer.parseInt(NminSetter.getText());
            int max = Integer.parseInt(NmaxSetter.getText());

            XYChart.Series err = tmp.TotalError(x0, y0, x1, min, max);
            TotalErrorChart.getData().addAll(err);
        }
    }

    private AtomicReference<Double> x0 = new AtomicReference<>((double) 1);
    private AtomicReference<Double> x1 = new AtomicReference<>((double) 6);
    private AtomicReference<Double> y0 = new AtomicReference<>((double) 1);
    private AtomicReference<Integer> N = new AtomicReference<>((200));

    void setValues(){
        XoSetter.setPromptText(String.valueOf(x0.get()));
        XoSetter.setText("");
        XSetter.setPromptText(String.valueOf(x1.get()));
        XSetter.setText("");
        YoSetter.setPromptText(String.valueOf(y0.get()));
        YoSetter.setText("");
        NSetter.setPromptText(String.valueOf(N.get()));
        NSetter.setText("");
    }

    void ClearGraphs(){
        MyChart.getData().clear();
        ErrorChart.getData().clear();
        TotalErrorChart.getData().clear();
    }

    @FXML
    void initialize() {
        getEulerData(x0.get(), y0.get(),x1.get(),N.get());

        SetBtn.setOnAction(event -> {
            double check = x0.get();
            if (XoSetter.getText() != null && !XoSetter.getText().trim().isEmpty()){ check = Double.parseDouble(XoSetter.getText()); }

            if(check != 0){
                if (XoSetter.getText() != null && !XoSetter.getText().trim().isEmpty()) x0.set(Double.parseDouble(XoSetter.getText()));
                if (XSetter.getText() != null && !XSetter.getText().trim().isEmpty()) x1.set(Double.parseDouble(XSetter.getText()));
                if (YoSetter.getText() != null && !YoSetter.getText().trim().isEmpty()) y0.set(Double.parseDouble(YoSetter.getText()));
                if (NSetter.getText() != null && !NSetter.getText().trim().isEmpty()) N.set(Integer.parseInt(NSetter.getText()));

                setValues();
                ClearGraphs();

                if (EulerChk.isSelected()) getEulerData(x0.get(), y0.get(),x1.get(),N.get());
                if (ImEulChk.isSelected()) getImprovedEulerData(x0.get(), y0.get(),x1.get(),N.get());
                if (RunKutChk.isSelected()) getRungeKuttaData(x0.get(), y0.get(),x1.get(),N.get());
            } else {
                XoSetter.setText("");
            }
        });
    }

    public void handleEulerChk(javafx.event.ActionEvent actionEvent) {
        EulerChk.setSelected(true);
        ImEulChk.setSelected(false);
        RunKutChk.setSelected(false);

        ClearGraphs();
        getEulerData(x0.get(), y0.get(),x1.get(),N.get());
        setValues();
    }

    public void handleImEulChk(ActionEvent event) {
        ImEulChk.setSelected(true);
        EulerChk.setSelected(false);
        RunKutChk.setSelected(false);

        ClearGraphs();
        getImprovedEulerData(x0.get(), y0.get(),x1.get(),N.get());
        setValues();
    }

    public void handleRunKutChk(ActionEvent event) {
        RunKutChk.setSelected(true);
        ImEulChk.setSelected(false);
        EulerChk.setSelected(false);

        ClearGraphs();
        getRungeKuttaData(x0.get(), y0.get(),x1.get(),N.get());
        setValues();
    }

    public void handleNChk(ActionEvent event) {

    }
}
