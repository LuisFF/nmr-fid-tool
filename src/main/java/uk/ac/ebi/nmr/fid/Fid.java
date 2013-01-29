package uk.ac.ebi.nmr.fid;

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
        data=new int[fid.size()];
        real=new int[fid.size()/2];
        imaginary=new int[fid.size()/2];

        for(int i=0; i < fid.size(); i++){
            data[i]=fid.get(i);
            if(i % 2 == 0)
                real[i/2]=data[i];
            else
                imaginary[(i-1)/2]=data[i];
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
