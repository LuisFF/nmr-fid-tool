package uk.ac.ebi.nmr.execs;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Fid;
import uk.ac.ebi.nmr.fid.FourierTransformTool;
import uk.ac.ebi.nmr.fid.io.AcquReader;
import uk.ac.ebi.nmr.fid.io.FidReader;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 10/01/2013
 * Time: 16:27
 * To change this template use File | Settings | File Templates.
 */
public class PlotFID02 {


        private XYSeriesCollection dataset;

        public static void main (String[] args) {
            new PlotFID02();
        }

        public PlotFID02() {

            dataset = new XYSeriesCollection();
            final XYSeries data = new XYSeries("data");
            Acqu acquisition = null;

            try {
                acquisition = new AcquReader("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                        "resources/examples/file_formats/bmse000109/1H/acqu").read();
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

            File fidFile = new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                    "resources/examples/file_formats/bmse000109/1H/fid");
            Fid fid = null;
            double[] spectrum = null;
            FourierTransformTool ft;
            try {
                fid = new FidReader(fidFile, acquisition).read();
                ft = new FourierTransformTool(fid,acquisition);

                spectrum = ft.computeFTT();
//                ft.fft(true);
//                spectrum=ft.getData();
                for (int i =0 ; i< fid.getData().length; i++){
                    data.add(i*ft.getProcessing().getDwellTime(), spectrum[i]);
                }
//            data.add(3, 2); //Point 1
//            data.add(1, 1); //Point 2
//            data.add(4, 1); //Point 3
//            data.add(2, 2); //Point 4
                dataset.addSeries(data);
                showGraph();
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }


        }

        public void showGraph() {
            final JFreeChart chart = createChart(dataset);
            final ChartPanel chartPanel = new ChartPanel(chart);
            chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
            final ApplicationFrame frame = new ApplicationFrame("Title");
            frame.setContentPane(chartPanel);
            frame.pack();
            frame.setVisible(true);
        }

        public JFreeChart createChart(final XYDataset dataset) {
            final JFreeChart chart = ChartFactory.createScatterPlot(
                    "Putrescine raw data from FID in time domain",                  // chart title
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

