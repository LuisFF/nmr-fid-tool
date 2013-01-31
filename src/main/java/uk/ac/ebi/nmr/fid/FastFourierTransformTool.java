/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.nmr.fid;

import org.apache.log4j.Logger;

/**
 * @name    FastFourierTransformTool
 * @date    2013.01.31
 * @version $Rev$ : Last Changed $Date$
 * @author  pmoreno
 * @author  $Author$ (this version)
 * @brief   ...class description...
 *
 */
public interface FastFourierTransformTool {

    double[] computeFFT();
    Proc getProcessing();
    
}
