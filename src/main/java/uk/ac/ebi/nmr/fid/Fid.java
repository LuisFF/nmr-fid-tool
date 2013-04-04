package uk.ac.ebi.nmr.fid;

import org.apache.commons.lang3.ArrayUtils;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:00
 * To change this template use File | Settings | File Templates.
 */
public class Fid {
    private int[] data;
    private int[] real;
    private int[] imaginary;

    public Fid(List<Integer> fid) {
        this(ArrayUtils.toPrimitive(fid.toArray(new Integer[fid.size()])));
    }

    public Fid(float [] fid) {
        int[]  fidInt = new int[fid.length];
        for(int i = 0 ; i< fid.length; i++){
            fidInt[i]=(int) Math.round(fid[i]);
        }
        new Fid(fidInt);

    }

    public Fid(int [] fid) {
        data = fid ;

        real=new int[fid.length/2];
        imaginary=new int[fid.length/2];

        for(int i=0; i < fid.length; i+=2){
            real[i/2]=data[i];// real are in even positions
            imaginary[i/2]=data[i+1];// imaginary are in odd positions
        }
    }

    public int[] getData() {
        return data;
    }

    public int[] getReal() {
        return real;
    }

    public int[] getImaginary() {
        return imaginary;
    }
}
