package uk.ac.ebi.nmr.fid.io;

import uk.ac.ebi.nmr.fid.Fid;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 21/03/2013
 * Time: 11:47
 * To change this template use File | Settings | File Templates.
 */
public interface FidReader {

    Fid read () throws Exception;

}
