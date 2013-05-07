package uk.ac.ebi.nmr.fid.tools.phasing;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
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
import uk.ac.ebi.nmr.fid.Spectrum;
import uk.ac.ebi.nmr.fid.io.BrukerAcquReader;
import uk.ac.ebi.nmr.fid.io.ConnjurFidReader;
import uk.ac.ebi.nmr.fid.io.FidReader;
import uk.ac.ebi.nmr.fid.tools.apodization.AbstractApodizator;
import uk.ac.ebi.nmr.fid.tools.apodization.ExponentialApodizator;

import java.io.File;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/05/2013
 * Time: 15:02
 * To change this template use File | Settings | File Templates.
 */
public class SimplePhaseCorrectorTest {

    XYSeries data;

    @Test
    public void testPhaseCorrection() throws Exception {
        data = new XYSeries("spectra");
        File fidFile = new File("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                "resources/examples/file_formats/bmse000109/1H/");
        Acqu acquisition = new BrukerAcquReader("/Users/ldpf/SVN/ldpf/dev/nmr-tools/src/test/java/"+
                "resources/examples/file_formats/bmse000109/1H/acqu").read();
        FidReader fidReader = new ConnjurFidReader(fidFile,acquisition);
        Spectrum spepSpectrum = fidReader.read();
        Assert.assertNotNull("fid was not properly read", spepSpectrum);
        double [] partialFID = new double[spepSpectrum.getFid().length-2*80];
        System.arraycopy(spepSpectrum.getFid(),160,partialFID,0,spepSpectrum.getFid().length-2*80);
        Spectrum partialSpectrum = new Spectrum(partialFID,spepSpectrum.getAcqu(),spepSpectrum.getProc());
        partialSpectrum.getAcqu().setAquiredPoints(spepSpectrum.getAcqu().getAquiredPoints()-2*80);
        DSPPhaseCorrection dspPhaseCorrection = new DSPPhaseCorrection();
        partialSpectrum= dspPhaseCorrection.dspPhaseCorrection(partialSpectrum);
        AbstractApodizator apodizator = new ExponentialApodizator(partialSpectrum);
        partialSpectrum = apodizator.calculate(0.1);

        DoubleFFT_1D fftd = new DoubleFFT_1D(partialSpectrum.getFid().length/2);
//            double[] spectrumRaw = apodizator.calculate(1.0);
        System.out.println(partialSpectrum.getFid().length);
        double [] realTransformed = partialSpectrum.getFid();
         double [] imgTransformed = partialSpectrum.getFid();
//        fftd.complexInverse(imgTransformed,true);
        fftd.realInverse(realTransformed,true);
        System.out.println(partialSpectrum.getFid().length);
//        fftd.complexInverse(spectrumRaw, true);
        double [] realChannel = new double[realTransformed.length/2];
        double [] imaginaryChannel = new double[realTransformed.length/2];
        // extract the real and imaginary side


        // extract the actual spectra from the quadrature spectrum
        for (int i = realTransformed.length/2; i< realTransformed.length; i+=2){
            realChannel[i/2]= realTransformed[2 * (i - realTransformed.length / 2)];
            imaginaryChannel[i/2]=realTransformed[2 * (i - realTransformed.length / 2)+1];
        }
        for (int i = 0; i< (realTransformed.length/2); i+=2){
            realChannel[i/2]= realTransformed[2*(i + realTransformed.length / 4)];
            imaginaryChannel[i/2]=realTransformed[2 * (i + realTransformed.length / 4)+1];
//            realChannel[i]= partialFID[2 * (i + partialFID.length / 2)];
//            imaginaryChannel[i/2]=partialFID[2 * (i + partialFID.length / 2)+1];
        }

        for(int i =0 ; i< realChannel.length; i++){
            data.add(i/2,realChannel[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection(data);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-spectra.png"), chart, 864, 1152);
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
