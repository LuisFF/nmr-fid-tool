package uk.ac.ebi.nmr.fid.io;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Fid;
import uk.ac.ebi.nmr.fid.Proc;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA. User: ldpf Date: 14/01/2013 Time: 14:12 To change this template use File | Settings |
 * File Templates.
 */
public class FidReader {

    private Acqu acquisition;
    private Proc processing;
    private InputStream fidInput;
    private double zeroFrequency = 0;

    public FidReader(InputStream fidInputS, Acqu acquisition) {
        this.acquisition = acquisition;
        this.fidInput = fidInputS;
    }

    public FidReader(InputStream fidInputS, Acqu acquisition, Proc processing) {
        this.acquisition = acquisition;
        this.fidInput = fidInputS;
        this.processing = processing;
    }

    public Fid read() throws Exception {

        Fid fid = null;

        DataInputStream in = null;

        in = new DataInputStream(new BufferedInputStream(fidInput));
        byte[] result = new byte[4];
        int[] fidInt = null;
        List<Integer> data = new ArrayList<Integer>();
        double[] fidDouble = null;
        int totalBytesRead = 0;
        System.out.println(acquisition.is32Bit());
        if (acquisition.is32Bit()) {
            fidInt = new int[acquisition.getAquiredPoints()];

            // read each set of 4 bytes in the fid to get the intensity values
            while (in.available() >= 4) {
                // byteOrder == 1 => Big-Endian else Little-Endian
                // this has to do with the order of the bytes
                data.add((acquisition.getByteOrder() == 1)
                        ? (in.readByte()) << 24
                        | (in.readByte() << 16)
                        | (in.readByte() << 8)
                        | (in.readByte())
                        : (in.readByte())
                        | (in.readByte() << 8)
                        | (in.readByte() << 16)
                        | (in.readByte() << 24));

            }

            System.out.println(data.isEmpty());
            System.out.println("bytes left: " + in.available());
            in.close();


            //split imaginary from real part
//                List<Integer>realPoints = new ArrayList<Integer>();
//                List<Integer> imagPoints = new ArrayList<Integer>();
//                for( int i = 0 ; i < points.size()/2 ; i++ ){
//                    realPoints.add(points.get(i*2)); //even index numbers contain real fid
//                    imagPoints.add(points.get(i*2+1));//odd index numbers contain imaginary fid
//                }

        } else {
            //TODO read 64bit FIDs... and store it as double?
            fidDouble = new double[acquisition.getAquiredPoints()];

        }

        System.out.println("Number of aquiredPoints and number of points in the FID:" + acquisition.getAquiredPoints() + " == "
                + data.size());

        System.out.println("Num bytes read: " + totalBytesRead);

        return new Fid(data);

    }
}
