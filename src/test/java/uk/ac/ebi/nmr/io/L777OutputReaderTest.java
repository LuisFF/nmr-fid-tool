package uk.ac.ebi.nmr.io;

import junit.framework.Assert;
import junit.framework.TestCase;
import org.junit.Test;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import java.util.HashMap;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 26/11/2012
 * Time: 11:29
 * To change this template use File | Settings | File Templates.
 */
public class L777OutputReaderTest extends TestCase {

    @Test
    public void testReadSerineVibrationalCorrections(){
        try {
            double [] chemShifts = {213.0158203, 115.1754230, 1.921899433,  -109.0486803,
                    110.2292075, 279.2940168, 108.2912881, 29.93661629, 30.43361457,
                    27.83464075, 27.32125633, 27.79092570, 28.95758297, 25.92665545};

            double [][] coupling = new double[14][14];
            coupling[7][8]=-7.68198;
            coupling[7][9]=5.81684;
            coupling[7][10]=1.87868;
            coupling[7][11]=-0.667782;
            coupling[7][12]=0.301501;
            coupling[7][13]=-0.264443;
            coupling[8][9]=12.6294;
            coupling[8][10]=-0.608863;
            coupling[8][11]=0.146698;
            coupling[8][12]= 0.0730429;
            coupling[8][13]=-0.268939;
            coupling[9][10]=4.43299;
            coupling[9][11]=9.75256;
            coupling[9][12]=-0.310240;
            coupling[9][13]=0.711682;
            coupling[10][11]=-10.0798;
            coupling[10][12]=12.9538;
            coupling[10][13]=-0.131227;
            coupling[11][12]=0.420182;
            coupling[11][13]=-0.161103;
            coupling[12][13]=-0.271349;


            L777OutputReader reader = new L777OutputReader(this.getClass()
                    .getResourceAsStream("/examples/file_formats/conf_10-vibcorr.log"));
            IChemFile chemFile = new ChemFile();
            chemFile = (IChemFile) reader.read(chemFile);
            // the reader will read a sequence of atomcontainer and only the last one will have all the information
            IAtomContainer atomContainer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                    ChemFileManipulator.getAllAtomContainers(chemFile).size()-1);

            for (int i = 0 ; i < atomContainer.getAtomCount();i++){
                IAtom atom = atomContainer.getAtom(i);
                Assert.assertEquals("Problem reading the magnetic tensor values ", chemShifts[i],
                        (Double) atom.getProperty(L777OutputReader.MAGNETIC_TENSOR));
            }

            // test the coupling between the hydrogens
            for (int i =7; i< atomContainer.getAtomCount();i++){
                Assert.assertNotNull(atomContainer.getAtom(i).getProperty(L777OutputReader.COUPLING_CONSTANTS));
                for (int j = i+1; j < atomContainer.getAtomCount(); j++){
                     HashMap<IAtom, Double> couplings =(HashMap<IAtom, Double>) atomContainer.getAtom(j)
                             .getProperty(L777OutputReader.COUPLING_CONSTANTS);
                    // get the coupling of atom j with the atom i
                    Assert.assertEquals("Problem reading the proton coupling constantes ",
                            coupling[i][j],couplings.get(atomContainer.getAtom(i)));
//                    System.out.println("Problem reading the proton coupling constantes "+
//                            coupling[i][j] +" " + couplings.get(atomContainer.getAtom(i)).getClass());
                }

            }

            //System.out.println(new SmilesGenerator().createSMILES(atomContainer).toString());
            String smile = "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]";
            Assert.assertEquals("AtomContainer has a different smile",
                    "[H]OC(=O)C([H])(N([H])[H])C([H])([H])O[H]",
                    new SmilesGenerator().createSMILES(atomContainer).toString());


        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }
}
