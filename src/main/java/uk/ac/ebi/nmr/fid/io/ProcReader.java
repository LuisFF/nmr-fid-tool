package uk.ac.ebi.nmr.fid.io;

import uk.ac.ebi.nmr.fid.Proc;

import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 02/04/2013
 * Time: 16:18
 * To change this template use File | Settings | File Templates.
 */
public interface ProcReader {

    Proc read() throws IOException;

}
