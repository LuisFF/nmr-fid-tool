package uk.ac.ebi.nmr.execs;

import edu.emory.mathcs.jtransforms.fft.FloatFFT_1D;
import edu.uchc.connjur.core.AnalysisDomain;
import edu.uchc.connjur.core.DimInfo;
import edu.uchc.connjur.spectrumtranslator.DataControl;
import edu.uchc.connjur.spectrumtranslator.SwDataSet;
import edu.uchc.connjur.spectrumtranslator.bruker.BrukerDataSetReader;
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
import uk.ac.ebi.nmr.fid.FastFourierTransformTool;
import uk.ac.ebi.nmr.fid.Fid;
import uk.ac.ebi.nmr.fid.JTransformFFTTool;
import uk.ac.ebi.nmr.fid.Proc;
import uk.ac.ebi.nmr.fid.io.AcquReader;
import uk.ac.ebi.nmr.fid.io.FidReader;
import uk.ac.ebi.nmr.fid.tools.PhaseCorrectionTool;

import java.io.File;
import java.io.InputStream;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA. User: ldpf Date: 10/01/2013 Time: 16:27 To change this template use File | Settings |
 * File Templates.
 */
public class PlotFID02 {

    private XYSeriesCollection dataset;

    public static void main(String[] args) throws Exception {
        new PlotFID02().showGraph();
    }

    public PlotFID02() throws Exception {

        dataset = new XYSeriesCollection();



//        initialApproach();
        connjurIntegration();


    }

    private void connjurIntegration() throws Exception {
        final XYSeries data = new XYSeries("data");
        BrukerDataSetReader dataSetReader = new BrukerDataSetReader();

//        dataSetReader.setDataSource(new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/" +
//                "resources/examples/file_formats/bmse000109/1H/"));
//        dataSetReader.setDataSource(new File("/Users/ldpf/Data/Databases/bmrb/tar_01/4_aminobenzoic_acid/nmr/bmse000066/1H/"));
        dataSetReader.setDataSource(new File("/Users/ldpf/Data/EIPOD/nmr/purecompound-201204-05/nmr/lf_gliotoxin_600/1/"));


        DataControl dataControl = new DataControl();
        SwDataSet dataSet = dataSetReader.read(null,dataControl);
//        SwDataSet dataSet = dataSetReader.read(true);

        System.out.println("Number of axes " + dataSet.getNumberAxes());

        System.out.println("Data type if the binary file: " + dataSet.getSourceDataType());
        DimInfo dimInfo = dataSet.getDataByDim(0);
        System.out.println("First order Phase :" + dimInfo.getFirstOrderPhase());
        System.out.println("Zero order Phase :" + dimInfo.getZeroOrderPhase());
        System.out.println("Total number of values (Re + Im):" + dimInfo.getTotalPoints());


        // get the frequency domain
        boolean isFrequency = false;
        for (int i = 0; i < dataSet.getNumberAxes() && !isFrequency; i++) {
            if (dataSet.getDataByDim(i).getAnalysisDomain() == AnalysisDomain.FREQUENCY) {
                isFrequency = true;
            }
        }
        // the real channel is accessible in the first fid?
        // the imaginary channel is accessible in the second fid?
        Iterator<edu.uchc.connjur.core.Fid> fidIterator = dataSet.fidIterator();
        if (fidIterator.hasNext()){
            // get the First fid - there should be more Fid for 2 dimensional spectra

            edu.uchc.connjur.core.Fid fidReal = fidIterator.next();
            //prepare for the appodization
            Acqu acquisition = new Acqu();
            acquisition.setAcquisitionMode(Acqu.AcquisitionMode.SIMULTANIOUS);
            System.out.println(acquisition.getAcquisitionMode().toString());
            acquisition.setAquiredPoints(fidReal.getValues().length);
            acquisition.setTransmiterFreq(500); // for bmse000109 is app. 500
            acquisition.setSpectralWidth(13); // for bmse000109 is app. 13
            Proc processing = new Proc(acquisition);
            processing.setLineBroadning(1);
            System.out.println("DW: " +processing.getDwellTime());
            System.out.println("LB: " +processing.getLineBroadning());


            int offSet = 141; // number of points ignored in the left side of the fid
            System.out.println("number of points in the re channel: " + (fidReal.getValues().length));
            edu.uchc.connjur.core.Fid fidIm = fidIterator.next();
            System.out.println("number of points in the im channel: "+ (fidIm.getValues().length));
            float [] fidTransformed= new float[(fidReal.getSize()-141)*2];
//            float [] fidImValues= new float[fidIm.getSize()-offSet];
            // interlace the real and imaginary points into a single array
            for (int i = 141; i < fidReal.getSize(); i++ ){
                fidTransformed[(i-141)*2]=fidReal.getValue(i);// insert the real (odd entries)
                fidTransformed[(i-141)*2+1]=fidIm.getValue(i);// insert the imaginary (even entries)
            }
            FloatFFT_1D fftd = new FloatFFT_1D(fidReal.getValues().length-141);

            fftd.complexInverse(fidTransformed, true);
            // reorganize the spectrum
            // swapt the left with the rigth side and invert the right side
            // do this at the global scale so that it can be used for the phase correction
            //TODO find the exact place to split the spectrum
//fftMirrored[1:((length(fftValue)+1)/2)]=fftValue[((length(fftValue)+1)/2):length(fftValue)];
//fftMirrored[((length(fftValue)+1)/2):length(fftValue)]=fftValue[1:((length(fftValue)+1)/2)];
            float [] spectrumRaw = new float[fidTransformed.length];

            for (int i = 0; i< fidTransformed.length/2; i++){ // check how to do with odd lengths...
                spectrumRaw[fidTransformed.length/2+i]=fidTransformed[i];
                spectrumRaw[i]=fidTransformed[fidTransformed.length/2+i];
            }


//            fftd.realInverse(fidImValues,true);
            // create the quadrature values for the spectra
//            quadLeft=Re(realPart)-Im(imagPart);
//            quadRight=Re(realPart)+Im(imagPart);






            // do the appodization
//            ApodizationTool apodizationTool = new ApodizationTool(convertFloatsToDoubles(fidRealValues),
//                    Acqu.AcquisitionMode.SIMULTANIOUS,processing);
//            double [] appodizedValues =apodizationTool.exponential();
//            double [] appodizedValues =convertFloatsToDoubles(fidRealValues);
            // prepare for the FFT
//            DoubleFFT_1D dfftd = new DoubleFFT_1D(appodizedValues.length/2);
//            dfftd.complexForward(appodizedValues);

            // do the phase correction
//            PhaseCorrectionTool phaseCorrectionTool = new PhaseCorrectionTool(appodizedValues, new Proc(acquisition));
//            double [] spectrumPhased = phaseCorrectionTool.zeroOrderPhasing(dimInfo.getZeroOrderPhase());
//            double [] spectrumPhased = phaseCorrectionTool.firstOrderPhasing(dimInfo.getFirstOrderPhase());
            // get just the real part of the points
//            double[] realPart = new double[appodizedValues.length/2];
//            int realIndex = 0;
//            for (int i = 0; i < appodizedValues.length-1; i+=2) {
//                realPart[realIndex++] = appodizedValues[i];
//            }
//            // put everything in data
//            for (int i=realPart.length-1; i>=0 ; i--){
//                data.add(realPart.length-i,realPart[i]);
//            }
            System.out.println("Total point of fid transformed (re+im)"+fidTransformed.length);
            // get only the real values, which should be in the odd entries
            //TODO break the spectra into half and mirror the left side...
            for (int i = 0; i< spectrumRaw.length/2; i++){
                data.add(i,spectrumRaw[i*2]);
            }

        }


//            data.add(3, 2); //Point 1
//            data.add(1, 1); //Point 2
//            data.add(4, 1); //Point 3
//            data.add(2, 2); //Point 4
        dataset.addSeries(data);
    }

    private void initialApproach() throws Exception {

        final XYSeries data = new XYSeries("data");
        Acqu acquisition = (new AcquReader(PlotFID02.class.getResourceAsStream("/examples/file_formats/bmse000109/1H/acqu"))).read();
        InputStream fidInput = PlotFID02.class.getResourceAsStream("/examples/file_formats/bmse000109/1H/fid");
        FastFourierTransformTool ft;

        Fid fid = new FidReader(fidInput, acquisition).read();
        System.out.println("numeber of points: " + fid.getData().length);

        ft = new JTransformFFTTool(fid, acquisition);
        System.out.println("total number of effective points: "+ ft.getProcessing().getTdEffective());

        Long milis = System.currentTimeMillis();
        double[] spectrum = ft.computeFFT();
        System.out.println("FFT timimg "+ft.getClass().getCanonicalName()+" "+(System.currentTimeMillis() - milis)+" ms ");
//                ft.fft(true);
//                spectrum=ft.getData();
        //for (int i = 0; i < fid.getData().length; i++) {

        PhaseCorrectionTool phaseCorrectionTool = new PhaseCorrectionTool(spectrum, new Proc(acquisition));
        double [] spectrumPhased = phaseCorrectionTool.zeroOrderPhasing(0.2);
        for (int i = 1; i < spectrum.length; i+=2) {
            //data.add(i * ft.getProcessing().getDwellTime(), spectrum[i]);
            data.add(i*ft.getProcessing().getDwellTime()/(acquisition.getAquiredPoints()*ft.getProcessing().getDwellTime()), spectrum[i]);
        }
        dataset.addSeries(data);
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
                "Gliotoxin 1H spectra", // chart title
                "Hz?", // x axis label
                "Intensity", // y axis label
                dataset, // data
                PlotOrientation.HORIZONTAL.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
                );
        XYPlot plot = (XYPlot) chart.getPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0, false);


        plot.setRenderer(renderer);
        return chart;
    }
    public static double[] convertFloatsToDoubles(float[] input)
    {
        if (input == null)
        {
            return null; // Or throw an exception - your choice
        }
        double[] output = new double[input.length];
        for (int i = 0; i < input.length; i++)
        {
            output[i] = input[i];
        }
        return output;
    }
}
