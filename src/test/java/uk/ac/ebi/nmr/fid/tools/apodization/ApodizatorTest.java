package uk.ac.ebi.nmr.fid.tools.apodization;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 22/01/2013
 * Time: 16:38
 * To change this template use File | Settings | File Templates.
 */
public class ApodizatorTest {

    static Acqu acquisition;
    static Proc processing;

    static double [] fid;
    static double [] fidObserved;
    static double [] fidApodizedEM;
    static double [] fidApodizedGM;
    static double [] fidApodizedLGM;
    static double [] fidApodizedTRAF;
    static double [] fidApodizedTRAFS;
    static double [] realFTObserved;
    static double [] imaginaryFTObserved;

    static double dw=0.01;

    XYSeries data;
    static XYSeries dataR;

    @BeforeClass
    public static void loadExternalData () {

        final Pattern DATAROW_PATTERN = Pattern.compile("(-?\\d+\\.?\\d*e?-?\\+?\\d*)");
        final Pattern DATA_PATTERN = Pattern
                .compile("(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)"+
                "\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)"+
                "\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)");
        
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(ApodizatorTest.class.getClassLoader()
                .getResourceAsStream("data/fid-simulated.csv")));
        int lines = 0;
        // count the number of lines
        // which corresponds to the number of data points if one removes the header line
        try {

            bufferedReader.readLine();
            while (bufferedReader.readLine() != null) lines++;

            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        fid = new double[lines];
        fidObserved = new double[lines];
        fidApodizedEM = new double[lines];
        fidApodizedGM = new double[lines];
        fidApodizedLGM = new double[lines];
        fidApodizedTRAF = new double[lines];
        fidApodizedTRAFS = new double[lines];
        realFTObserved = new double[lines];
        imaginaryFTObserved = new double[lines];

        bufferedReader = new BufferedReader(new InputStreamReader(ApodizatorTest.class.getClassLoader()
                .getResourceAsStream("data/fid-simulated.csv")));



        try {
            String line = bufferedReader.readLine();
            int index = 0;
            Matcher matcher = null;
            dataR = new XYSeries("dataR");

            while (line != null){
                matcher = DATAROW_PATTERN.matcher(line);
                if(matcher.find()){
                    int matchEnd=0;
                    int times =0;
                    while (matcher.find(matchEnd)){
                        matchEnd=matcher.end();
//                        System.out.println("Pattern found "+ new String(new char[times]).replace("\0","\t" )
//                                +Double.parseDouble(matcher.group(1)));
                        switch (times){
                            case 0 :
                                fid[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 1 :
                                fidObserved[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 2 :
                                realFTObserved[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 3 :
                                imaginaryFTObserved[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 4 :
                                fidApodizedEM[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 5 :
                                fidApodizedGM[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 6 :
                                fidApodizedLGM[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 7 :
                                fidApodizedTRAF[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 8 :
                                fidApodizedTRAFS[index]=Double.parseDouble(matcher.group(1));
                                break;
                            default:
                                break;
                        }
                        ++times;
                    }
                    dataR.add(index*dw,fidObserved[index]);
                    ++index;
                }

                line=bufferedReader.readLine();
            }
            System.out.println(index);
            acquisition = new Acqu();
            // The test data I have is sequencial, meaning that I have only the FID of the real channel
            // The simultaneous would have the value of the real channel followed by the corresponding value of the
            // imaginary chanel
            acquisition.setAcquisitionMode(Acqu.AcquisitionMode.SEQUENTIAL);
            System.out.println(acquisition.getAcquisitionMode().toString());
            acquisition.setAquiredPoints(fid.length);
            acquisition.setTransmiterFreq(8 + 1 / 3); // for bmse000109 is app. 500
            acquisition.setSpectralWidth(6.25); // for bmse000109 is app. 13
            processing = new Proc(acquisition);
            processing.setLineBroadening(1);
            processing.setGbFactor(0.03);
            System.out.println("DW: " +processing.getDwellTime());
            System.out.println("LB: " +processing.getLineBroadening());
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        }

        @Test
        public void TestExponentialApodizator(){
            data = new XYSeries("data");
            Apodizator expApodizator= new ExponentialApodizator(fidObserved,
                    Acqu.AcquisitionMode.SEQUENTIAL,processing);

            double[] apodizedExponential= new double[0];
            try {
                apodizedExponential = expApodizator.calculate(1);
            boolean controlSTERR = false;
            for (int i =0; i< fidObserved.length;i++){
                data.add(i*dw,apodizedExponential[i]);
                Assert.assertTrue("Numerical deviation above 1E-12",
                        Math.abs(apodizedExponential[i]-fidApodizedEM[i])<1E-12);

                if(!(controlSTERR) && Math.abs(apodizedExponential[i]-fidApodizedEM[i])>1E-12) {
                    System.err.println("Exponential apodization with deviation larger than 1E-12");
                    System.err.println(i+" "+ fidApodizedEM[i]+" "+ apodizedExponential[i]);
                    controlSTERR=true;
                }
            }
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            XYSeriesCollection dataset = new XYSeriesCollection(dataR);
            JFreeChart chart = createChart(dataset);
            dataset = new XYSeriesCollection(data);
            JFreeChart chart2 = createChart(dataset);
            try {
                ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart.png"), chart, 864 ,1152);
                ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart2.png"), chart2, 864 ,1152);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
    }


    @Test
    public void TestGaussianApodizator(){
        data = new XYSeries("data");
        Apodizator gaussApodizator= new GaussianApodizator(fidObserved,
                Acqu.AcquisitionMode.SEQUENTIAL,processing);
        double[] apodized = new double[0];
        try {
            apodized = gaussApodizator.calculate(1);
            boolean controlSTERR = false;
            for (int i =0; i< fidObserved.length;i++){
                data.add(i*dw,apodized[i]);
                Assert.assertTrue("Numerical deviation above 1E-12",
                        Math.abs(apodized[i]-fidApodizedGM[i])<1E-12);
                if(Math.abs(apodized[i]-fidApodizedGM[i])>1E-12) {
                    if(!controlSTERR){
                        System.err.println("Gaussian apodization with deviation larger than 1E-12");
                        System.err.println(i+" "+ fidApodizedGM[i]+" "+ apodized[i]);
                        controlSTERR=true;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        JFreeChart chart = createChart(dataset);
        dataset = new XYSeriesCollection(data);
        JFreeChart chart2 = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart.png"), chart, 864 ,1152);
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart2.png"), chart2, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
    @Test
    public void TestLorentzGausApodizator(){
        data = new XYSeries("data");
        Apodizator lgApodizator= new LorentzGausApodizatior(fidObserved,
                Acqu.AcquisitionMode.SEQUENTIAL,processing);

        double[] apodized= new double[0];
        try {
            apodized = lgApodizator.calculate(-0.3);
            boolean controlSTERR = false;
            int count=0;
            for (int i =0; i< fidObserved.length;i++){
                data.add(i*dw,apodized[i]);

                if(Math.abs(apodized[i]-fidApodizedLGM[i])>1E-3) {
                    if(!(controlSTERR)){
                    System.err.println("Lorentz-Gauss apodization with deviation larger than 1E-3");
                    System.err.println(i+" "+fidApodizedLGM[i]+" "+ apodized[i]);
                        controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("Lorentz-Gauss apodization deviations above 1E-3: "+count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        JFreeChart chart = createChart(dataset);
        dataset = new XYSeriesCollection(data);
        JFreeChart chart2 = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart.png"), chart, 864 ,1152);
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart2.png"), chart2, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Test
    public void TestTRAFApodizator(){
        data = new XYSeries("data");
        Apodizator tarfApodizator= new TrafApodizator(fidObserved,
                Acqu.AcquisitionMode.SEQUENTIAL,processing);

        double[] apodized= new double[0];
        try {
            apodized = tarfApodizator.calculate(0.06);
            boolean controlSTERR = false;
            int count = 0;
            for (int i =0; i< fidObserved.length;i++){
                data.add(i*dw,apodized[i]);
                Assert.assertTrue("Numerical deviation above 1E-3",
                        Math.abs(apodized[i]-fidApodizedTRAF[i])<1E-3);
                if(Math.abs(apodized[i]-fidApodizedTRAF[i])>1E-4) {
                    if(!(controlSTERR)){
                    System.err.println("TRAF apodization with deviation larger than 1E-4");
                    System.err.println(i+" "+fidApodizedTRAF[i]+" "+ apodized[i]);
                    controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("TRAF apodization deviations above 1E-4: "+count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        JFreeChart chart = createChart(dataset);
        dataset = new XYSeriesCollection(data);
        JFreeChart chart2 = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart.png"), chart, 864 ,1152);
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart2.png"), chart2, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Test
    public void TestTRAFSApodizator(){
        data = new XYSeries("data");
        Apodizator tarfsApodizator= new TrafsApodizator(fidObserved,
                Acqu.AcquisitionMode.SEQUENTIAL,processing);

        double[] apodized= new double[0];
        try {
            apodized = tarfsApodizator.calculate(0.06);
            boolean controlSTERR = false;
            int count = 0;
            for (int i =0; i< fidObserved.length;i++){
                data.add(i*dw,apodized[i]);
                Assert.assertTrue("Numerical deviation above 1E-3",
                        Math.abs(apodized[i]-fidApodizedTRAFS[i])<1E-3);
                if(Math.abs(apodized[i]-fidApodizedTRAFS[i])>1E-4) {
                    if(!(controlSTERR) ){
                    System.err.println("TRAFS apodization with deviation larger than 1E-4");
                    System.err.println(i+" "+fidApodizedTRAFS[i]+" "+ apodized[i]);
                    controlSTERR=true;
                    }
                    count++;
                }
            }
            System.out.println("TRAFS apodization deviations above 1E-4: "+count);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        XYSeriesCollection dataset = new XYSeriesCollection(dataR);
        JFreeChart chart = createChart(dataset);
        dataset = new XYSeriesCollection(data);
        JFreeChart chart2 = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart.png"), chart, 864 ,1152);
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart2.png"), chart2, 864 ,1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private JFreeChart createChart(final XYDataset dataset) {
        final JFreeChart chart = ChartFactory.createScatterPlot(
                "Fid",                  // chart title
                "Hz?",                      // x axis label
                "Intensity",                      // y axis label
                dataset,                  // data
                PlotOrientation.HORIZONTAL.VERTICAL,
                true,                     // include legend
                true,                     // tooltips
                false                     // urls
        );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0,false);


        plot.setRenderer(renderer);
        return chart;
    }

}
