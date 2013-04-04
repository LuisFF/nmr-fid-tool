package uk.ac.ebi.nmr.fid.tools.apodization;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 11:39
 * To change this template use File | Settings | File Templates.
 */
public interface Apodizator {

    double[] calculate() throws Exception;

    double[] calculate(double lineBroadning) throws Exception;
}
