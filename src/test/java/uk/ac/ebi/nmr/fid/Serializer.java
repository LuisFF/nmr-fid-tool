/*
 * Copyright (c) 2013 EMBL, European Bioinformatics Institute.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package uk.ac.ebi.nmr.fid;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import uk.ac.ebi.nmr.fid.io.FidReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to serialize the data obtained from R
 *
 * @author  Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 07/05/2013
 * Time: 17:15
 * To change this template use File | Settings | File Templates.
 */
public class Serializer {

    public static void serialize() throws IOException {
        /*
             serialize the fid data from R
              */

        // pattern for scientific notation
//        final Pattern DATAROW_PATTERN = Pattern.compile("\"(-?\\s?\\d+\\.\\d+e[-,\\+]\\d+)\"");
//        final Pattern DATAROW_PATTERN = Pattern.compile("(TRUE|FALSE)");
        // pattern for FIDs...
        final Pattern DATAROW_PATTERN = Pattern.compile("(-?\\d+\\.?\\d*)");
        BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(Serializer.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-fft-phased-1i.csv")));
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
//        boolean[] fid = new boolean[lines];
        bufferedReader = new BufferedReader(new InputStreamReader(FidReader.class.getClassLoader()
                .getResourceAsStream("data/bmse000109/1h/fid-fft-phased-1i.csv")));

        String line = bufferedReader.readLine();
        int index = 0;
        Matcher matcher = null;
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
//                            fid[index]=Boolean.parseBoolean(matcher.group(1));
                            data.add(index,fid[index]);
//                            data.add(index,fid[index]?1:0);
                            break;
                        default:
                            break;}
                    ++times;
                }

                ++index;
            }
            line=bufferedReader.readLine();
        }
        System.out.println(System.getProperty("os.name"));
    FileOutputStream fout = new FileOutputStream(new File("/Users/ldpf/SVN/ldpf/dev/nmr-fid-tool/src/test/java/" +
            "resources/data/bmse000109/1h/fid-fft-phased-1i.ser"));
    ObjectOutputStream oos = new ObjectOutputStream(fout);
    oos.writeObject(fid);
    oos.flush();
    oos.close();
    fout.close();
        XYSeriesCollection dataset = new XYSeriesCollection(data);
        JFreeChart chart = createChart(dataset);
        try {
            ChartUtilities.saveChartAsPNG(new File("/Users/ldpf/Downloads/chart-serialization.png"), chart, 864, 1152);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    public static void main(String[] args) {
        try {
            serialize();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private static JFreeChart createChart(final XYDataset dataset) {
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

    public static Object loadSerializedObject(InputStream inputStream) throws IOException, ClassNotFoundException {
        long startTime = System.currentTimeMillis();
        ObjectInputStream ois = new ObjectInputStream(inputStream);
        Object serializedOject = ois.readObject();
        long endTime   = System.currentTimeMillis();
        System.out.println("Time reading object in "+inputStream.toString()+":\t\t"+(endTime - startTime)+" ms");
        return serializedOject;
    }
}
