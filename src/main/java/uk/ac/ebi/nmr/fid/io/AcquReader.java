package uk.ac.ebi.nmr.fid.io;

import uk.ac.ebi.nmr.fid.Acqu;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:24
 * To change this template use File | Settings | File Templates.
 */
public interface AcquReader {
    Acqu read() throws Exception;
}
