package uk.ac.ebi.nmr.fid.tools;

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
public class ApodizationToolTest {

    @Test
    public void TestApodizationTool (){


        final Pattern DATAROW_PATTERN = Pattern.compile("(-?\\d+\\.?\\d*e?-?\\+?\\d*)");
        final Pattern DATA_PATTERN = Pattern
                .compile("(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)"+
                "\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)"+
                "\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)\\s+(-?\\d+\\.\\d+e?-?\\+?\\d+)");
        
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(this.getClass().getClassLoader()
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



        double [] fid = new double[lines];
        double [] fidObserved = new double[lines];
        double [] fidApodizedEM = new double[lines];
        double [] fidApodizedGM = new double[lines];
        double [] fidApodizedDE = new double[lines];
        double [] fidApodizedTRAF = new double[lines];
        double [] realFTObserved = new double[lines];
        double [] imaginaryFTObserved = new double[lines];


        bufferedReader = new BufferedReader(new InputStreamReader(this.getClass().getClassLoader()
                .getResourceAsStream("data/fid-simulated.csv")));

        double dw=0.01;

        try {
            String line = bufferedReader.readLine();
            int index = 0;
            Matcher matcher = null;
            XYSeries dataR = new XYSeries("dataR");
            XYSeries data = new XYSeries("data");
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
                                fidApodizedDE[index]=Double.parseDouble(matcher.group(1));
                                break;
                            case 7 :
                                fidApodizedTRAF[index]=Double.parseDouble(matcher.group(1));
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
            Acqu acquisition = new Acqu();
            acquisition.setAcquisitionMode(Acqu.AcquisitionMode.SIMULTANIOUS);
            System.out.println(acquisition.getAcquisitionMode().toString());
            acquisition.setAquiredPoints(fid.length);
            acquisition.setTransmiterFreq(8+1/3); // for bmse000109 is app. 500
            acquisition.setSpectralWidth(6.25); // for bmse000109 is app. 13
            Proc processing = new Proc(acquisition);
            processing.setLineBroadning(1);
            System.out.println("DW: " +processing.getDwellTime());
            System.out.println("LB: " +processing.getLineBroadning());

            ApodizationTool apodizationTool = new ApodizationTool(fidObserved,
                    Acqu.AcquisitionMode.SIMULTANIOUS,processing);



            double[] apodizedExponential= apodizationTool.exponential();
            double[] apodizedGaussian= apodizationTool.gaussianBramer2001(1);
            int count=0;
            boolean exponetialDeviated=false;
            boolean gaussianBramerDeviated=false;
            for (int i =0; i< fidObserved.length;i++){
                data.add(i*dw,apodizedExponential[i]);
                //check exponential apodization results
                Assert.assertTrue("Exponential apodization with a deviation larger than 0.01",
                       (Math.abs(apodizedExponential[i] - fidApodizedEM[i]) < 0.01));
                Assert.assertTrue("Gaussian apodization with a deviation larger than 0.01",
                        (Math.abs(apodizedGaussian[i] - fidApodizedGM[i]) < 0.01));

                if((! exponetialDeviated) && Math.abs(apodizedExponential[i]-fidApodizedEM[i])>0.001) {
                    System.err.println("Exponential apodization with deviation larger than 0.001");
                    System.err.println(fidApodizedEM[i]+" "+ apodizedExponential[i]);
                    exponetialDeviated=true;
                }
                if((! gaussianBramerDeviated) && Math.abs(apodizedGaussian[i]-fidApodizedGM[i])>0.001) {
                    System.err.println("Gaussian Brammer 2001 apodization with deviation larger than 0.001");
                    System.err.println(fidApodizedGM[i]+" "+ apodizedGaussian[i]);
                    gaussianBramerDeviated=true;
                }
//                if(Math.abs(apodizedGaussian[i]-fidApodizedGM[i])>0.001)
//                    System.out.println(fidApodizedGM[i]+" "+ apodizedGaussian[i]);
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


//
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
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
