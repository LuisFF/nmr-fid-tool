package uk.ac.ebi.nmr.io;

/* $Revision$ $Author$ $Date$
*
* Copyright (C) 2002-2003  Bradley A. Smith <yeldar@home.com>
*
* Contact: cdk-devel@lists.sourceforge.net
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*  This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
*/

import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.rebond.RebondTool;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.formats.Gaussian03Format;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.SaturationChecker;

import javax.vecmath.Point3d;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * A reader for Gaussian09 output.
 * Gaussian 09 is a quantum chemistry program
 * by Gaussian, Inc. (<a href="http://www.gaussian.com/">http://www.gaussian.com/</a>).
 * <p/>
 * <p>Molecular coordinates, energies, and normal coordinates of
 * vibrations are read. Each set of coordinates is added to the
 * ChemFile in the order they are found. Energies and vibrations
 * are associated with the previously read set of coordinates.
 * <p/>
 * <p>This reader was developed from a small set of
 * example output files, and therefore, is not guaranteed to
 * properly read all Gaussian09 output. If you have problems,
 * please contact the author of this code, not the developers
 * of Gaussian09.
 * <p/>
 * <p>This code was adaptated by Luis from Gaussian03Reader written by
 * Jonathan and Bradley, and ported to CDK by  Egon.
 * <p/>
 * Date: 13/08/2012
 *
 * @author Luis F. de Figueiredo
 * @author Jonathan C. Rienstra-Kiracofe <jrienst@emory.edu>
 * @author Bradley A. Smith <yeldar@home.com>
 * @author Egon Willighagen
 * @cdk.module io
 * @cdk.githash
 */

public class G09OutputReader extends DefaultChemObjectReader {

    public static final String BOND_LENGTH = "org.openscience.cdk.io.G09OutputReader:bond.length";
    // energy is in hartrees: Eh = 4.359 744 34 x 10-18 J
    public static final String STRUCTURE_ENERGY = "org.openscience.cdk.io.G09OutputReader:structure.energy";
    public static final String OPTIMIZATION_CYCLES = "org.openscience.cdk.io.G09OutputReader:optimization.cycles";
    public static final String MAGNETIC_TENSOR = "org.openscience.cdk.io.G09OutputReader:isotropic.magnetic.tensor";
    public static final String COUPLING_CONSTANTS = "org.openscience.cdk.io.G09OutputReader:isotropic.coupling.constants";
    public static final String NORMAL_TERMINATION = "org.openscience.cdk.io.G09OutputReader:normal.termination";

    // Tables required to define the IAtomContainer in the input file
    private final static Pattern REGEX_Z_MATRIX = Pattern.compile("Symbolic\\sZ-matrix:");
    private final static Pattern REGEX_CONNECTION_TABLE = Pattern.compile("Initial Parameters");
    // assume that at least there is one digit in the left of the decimal place
    private final static Pattern REGEX_INIT_SET_COORDINATES = Pattern
            .compile("\\s?(\\w+|\\d+)\\s+(-?\\d+\\.\\d*)\\s+(-?\\d+\\.\\d*)\\s+(-?\\d+\\.\\d*)");
    private final static Pattern REGEX_BOND_EDGE = Pattern
            .compile("R\\d+\\s+R\\((\\d+),(\\d+)\\)\\s+?(-?\\d+\\.\\d*)");

    private final static Pattern REGEX_END_TABLE= Pattern.compile("------------------");

    // geometry table and parameters
    private final static Pattern REGEX_GEOMETRY_START = Pattern.compile("Input orientation:");
    private final static Pattern REGEX_OPTIMIZED_GEOMETRY = Pattern.compile("\\s+(\\d+)\\s+(\\d+)\\s+\\d+\\s+"+
                            "(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)");

    private final static Pattern REGEX_SCF_ENERGY = Pattern.compile("SCF Done:.+=\\s*(-?\\d+\\.\\d*).+(\\d+)\\scycles");

    // NMR tables
    private final static Pattern REGEX_NMR_PARAMETERS_START = Pattern.compile("(Calculating GIAO)|(NMR Shielding)");
    protected final static Pattern REGEX_NMR_PARAMETERS_STOP = Pattern.compile("(\\*\\*\\*\\*\\*)|(==========)");

    protected final static Pattern REGEX_SS_COUPLING_FLAG = Pattern.compile("Total nuclear spin-spin coupling J");
    private final static Pattern REGEX_MAGNETIC_TENSOR_TABLE = Pattern.compile(
            "\\s+(\\d+)\\s+(\\w+)\\s+Isotropic\\s=\\s+(-?\\d+\\.\\d*(D-?\\+?\\d{2})?)");
//    private final static Pattern REGEX_MAGNETIC_TENSOR_TABLE02 = Pattern.compile(
//            "\\s+(\\d+)\\s+(\\w+)\\s+Isotropic\\s=\\s+(-?\\d+\\.\\d{10}D-?\\+?\\d{2})");

    protected final static Pattern REGEX_SS_COUPLING_HEADER_FLAG = Pattern
            .compile("\\s{8}");
    protected final static Pattern REGEX_SS_COUPLING_DIAGONAL = Pattern
            .compile("\\s+(\\d+)\\s+(-?\\d+\\.\\d{6}D-?\\+?\\d{2})");
    protected final static Pattern REGEX_SS_COUPLING = Pattern
            .compile("\\s+(-?\\d+\\.\\d{6}D-?\\+?\\d{2})");

    private final static Pattern REGEX_LEVEL_OF_THEORY = Pattern.compile("GINC");

    private final static Pattern REGEX_NORMAL_TERMINATION = Pattern.compile("Normal termination of Gaussian");

    private static String levelOfTheory;

    protected BufferedReader input;
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(G09OutputReader.class);

    public G09OutputReader(Reader reader) {
        input = new BufferedReader(reader);
    }

    public G09OutputReader(InputStream input) {
        this(new InputStreamReader(input));
    }

    public G09OutputReader() {
        this(new StringReader(""));
    }

    @TestMethod("testGetFormat")
    public IResourceFormat getFormat() {
        return Gaussian03Format.getInstance();
    }

    @TestMethod("testSetReader_Reader")
    public void setReader(Reader reader) throws CDKException {
        this.input = new BufferedReader(input);
    }

    @TestMethod("testSetReader_InputStream")
    public void setReader(InputStream input) throws CDKException {
        setReader(new InputStreamReader(input));
    }

    @TestMethod("testAccepts")
    public boolean accepts(Class classObject) {
        Class[] interfaces = classObject.getInterfaces();
        for (int i = 0; i < interfaces.length; i++) {
            if (IChemFile.class.equals(interfaces[i])) return true;
            if (IChemSequence.class.equals(interfaces[i])) return true;
        }
        return false;
    }

    public <T extends IChemObject> T read(T object) throws CDKException {
        if (object instanceof IChemSequence) {
            return (T) readChemSequence((IChemSequence) object);
        } else if (object instanceof IChemFile) {
            return (T) readChemFile((IChemFile) object);
        } else {
            throw new CDKException("Object " + object.getClass().getName() + " is not supported");
        }
    }

    @TestMethod("testClose")
    public void close() throws IOException {
        input.close();
    }

    private IChemFile readChemFile(IChemFile chemFile) throws CDKException {
        IChemSequence sequence = readChemSequence(chemFile.getBuilder().newInstance(IChemSequence.class));
        chemFile.addChemSequence(sequence);
        return chemFile;
    }

    /**
     * @param sequence
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    private IChemSequence readChemSequence(IChemSequence sequence) throws CDKException {

        IChemModel model = null;
        Set<IBond> bondsSet = null;
        Matcher matcher;


        try {
            String line = input.readLine();
            // Define the atomcontainer given as input
            // this is usually defined in the top of the file before the optimization start
            // the connection between the atoms requires though an initial optiomization
            // which should be done by gaussian (as default?)
            model = sequence.getBuilder().newInstance(IChemModel.class);
            // by default normal termination is false
            model.setProperty(NORMAL_TERMINATION,false);
            while (input.ready() && (line != null)) {
                if (REGEX_Z_MATRIX.matcher(line).find()) {
                    // Found a set of coordinates
                    try {
                        readAtomContainer(model);
                    } catch (IOException exception) {
                        throw new CDKException("Error while reading coordinates: " + exception.toString(), exception);
                    }
                    break;
                }
                line = input.readLine();
            }
            // extract additional information from the file
            if (model != null) {
                // Read all other data
                String lastLine = "";
                while (input.ready() && (line != null)) {
                    // update the geometry
                    if (REGEX_GEOMETRY_START.matcher(line).find()) {
                        // Found a new set of coordinates
                        // Add current frame to file and create a new one.
                        try {
                            sequence.addChemModel((IChemModel) model.clone());
                            fireFrameRead();
                            // carry on with the same model just update the coordinates
                            // and possible add more information
                            updateCoordinates(model);
                        } catch (CloneNotSupportedException e) {
                            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        } catch (Exception e) {
                            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                    }

                    if (REGEX_SCF_ENERGY.matcher(line).find()) {
                        // Found an energy
                        matcher = REGEX_SCF_ENERGY.matcher(line);
                        matcher.find();
                        model.setProperty(STRUCTURE_ENERGY, matcher.group(1));
                        model.setProperty(OPTIMIZATION_CYCLES, matcher.group(2));
                    }
                    if (line.indexOf("Harmonic frequencies") >= 0) {
                        // Found a set of vibrations
//                        try {
//                            readFrequencies(model);
//                        } catch (IOException exception) {
//                            throw new CDKException("Error while reading frequencies: " + exception.toString(), exception);
//                        }
                    }
                    if (line.indexOf("Mulliken atomic charges") >= 0) {
//                            readPartialCharges(model);
                    }
                    if (REGEX_NMR_PARAMETERS_START.matcher(line).find()) {
                        // Found NMR data
                        try {
                            readNMRData(model, line);
                        } catch (IOException exception) {
                            throw new CDKException("Error while reading NMR data: " + exception.toString(), exception);
                        } catch (Exception e) {
                            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                    }
                    if (REGEX_LEVEL_OF_THEORY.matcher(line).find()) {
                        // Found calculation level of theory
                        levelOfTheory = parseLevelOfTheory(line);
                    }
                    line = input.readLine();
                    lastLine=line;
                }
                // check if gaussian terminated correctly
                if (REGEX_NORMAL_TERMINATION.matcher(lastLine).find())
                    model.setProperty(G09OutputReader.NORMAL_TERMINATION, true);


                // Add current frame to file
                sequence.addChemModel(model);
                fireFrameRead();
            }
        } catch (IOException exception) {
            throw new CDKException("Error while reading general structure: " + exception.toString(), exception);
        }
        return sequence;
    }

    /**
     * read the file and extract the atomcontainer
     *
     * @param model
     */
    private void readAtomContainer(IChemModel model) throws CDKException, IOException {
        IAtomContainer container = model.getBuilder().newInstance(IAtomContainer.class);
        String line;
        Matcher matcher;
        boolean reachBound = false;
        int countEndTablePattern=0;
        int bondNb = 1;

        // default values suggested by Eagon
        RebondTool rebonder = new RebondTool(2.0,0.5,0.5);


        while (input.ready()) {
            line = input.readLine();
            if (REGEX_INIT_SET_COORDINATES.matcher(line).find()) {
                //in previous versions of the reader, the symbol was defined based in the atomicNumber
                matcher = REGEX_INIT_SET_COORDINATES.matcher(line);
                matcher.find();

                IAtom atom = model.getBuilder().newInstance(IAtom.class, matcher.group(1));

                atom.setPoint3d(new Point3d(Double.parseDouble(matcher.group(2)),
                        Double.parseDouble(matcher.group(3)),
                        Double.parseDouble(matcher.group(4))));
                container.addAtom(atom);
            }
            // read bonds from the first optimization
//            if (REGEX_BOND_EDGE.matcher(line).find()) {
//                matcher = REGEX_BOND_EDGE.matcher(line);
//                matcher.find();
//
//                IBond bond = new Bond(container.getAtom(Integer.parseInt(matcher.group(1)) - 1),
//                        container.getAtom(Integer.parseInt(matcher.group(2)) - 1),
//                        IBond.Order.SINGLE);
//                bond.setProperty(BOND_LENGTH, Double.parseDouble(matcher.group(3)));
//                container.addBond(bond);
//            }
            countEndTablePattern+=(REGEX_END_TABLE.matcher(line).find())?1:0;
            // this is a weak way of breaking the reading
            if ((line == null) || (countEndTablePattern > 1)) {
                break;
            }
        }

        // perform atom typing with AtomTypeFactory instead of AtomContainerManipulator
        // the latter doesn't work
        AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/config/data/jmol_atomtypes.txt"
                , container.getBuilder());
        for (IAtom atom : container.atoms())
                factory.configure(atom);

        rebonder.rebond(container);

        CDKHueckelAromaticityDetector.detectAromaticity(container);
        new SaturationChecker().saturate(container);

        // define the model as a molecule set
        // in order to store the different geometries
        IMoleculeSet moleculeSet = model.getBuilder().newInstance(IMoleculeSet.class);
        moleculeSet.addMolecule(model.getBuilder().newInstance(IMolecule.class, container));
        model.setMoleculeSet(moleculeSet);
    }

    /**
     * Updates the coordinates of an existing chemical model
     *
     * @param model the destination ChemModel
     * @throws java.io.IOException if an I/O error occurs
     */
    private void updateCoordinates(IChemModel model) throws Exception, IOException, CloneNotSupportedException {
        //TODO confirm if in the previous version one is creating a new model for each optimized geometry
        // or for each optimization procedure
        // clone the IAtomContainer in order to keep the existing properties
        IAtomContainer container = (IAtomContainer) model.getMoleculeSet()
                .getMolecule(model.getMoleculeSet().getMoleculeCount()-1).clone();
        String line = input.readLine();
        int countEndTablePattern=0;
        Matcher matcher;
        while (input.ready() & line !=null) {
            line = input.readLine();
            countEndTablePattern+=(REGEX_END_TABLE.matcher(line).find())?1:0;
            if ((line == null) || (countEndTablePattern > 1)){
                break;
            }
            if (REGEX_OPTIMIZED_GEOMETRY.matcher(line).find()){
                matcher=REGEX_OPTIMIZED_GEOMETRY.matcher(line);
                matcher.find();
                // simple test to check if at least atoms of with the same atomic number are in the same position
                if(!container.getAtom(Integer.parseInt(matcher.group(1))-1).getAtomicNumber()
                        .equals(Integer.parseInt(matcher.group(2))))
                    throw new Exception("The atoms are scrambled");
                // update the coordinates of the atom
                container.getAtom(Integer.parseInt(matcher.group(1))-1).setPoint3d(new Point3d(
                        Double.parseDouble(matcher.group(3)),
                        Double.parseDouble(matcher.group(4)),
                        Double.parseDouble(matcher.group(5))));
            }
        }

        IMoleculeSet moleculeSet = model.getBuilder().newInstance(IMoleculeSet.class);
        moleculeSet.addMolecule(model.getBuilder().newInstance(IMolecule.class, container));
        model.setMoleculeSet(moleculeSet);
    }

    /**
     * Reads partial atomic charges and add the to the given ChemModel.
     */
    private void readPartialCharges(IChemModel model) throws CDKException, IOException {
        //TODO convert this method to use regular expressions
        logger.info("Reading partial atomic charges");
        IAtomContainer molecule = (IAtomContainer) model.getMoleculeSet()
                .getMolecule(model.getMoleculeSet().getMoleculeCount()-1);
        String line = input.readLine(); // skip first line after "Total atomic charges"
        while (input.ready()) {
            line = input.readLine();
            logger.debug("Read charge block line: " + line);
            if ((line == null) || (line.indexOf("Sum of Mulliken charges") >= 0)) {
                logger.debug("End of charge block found");
                break;
            }
            StringReader sr = new StringReader(line);
            StreamTokenizer tokenizer = new StreamTokenizer(sr);
            if (tokenizer.nextToken() == StreamTokenizer.TT_NUMBER) {
                int atomCounter = (int) tokenizer.nval;

                tokenizer.nextToken(); // ignore the symbol

                double charge = 0.0;
                if (tokenizer.nextToken() == StreamTokenizer.TT_NUMBER) {
                    charge = (double) tokenizer.nval;
                    logger.debug("Found charge for atom " + atomCounter +
                            ": " + charge);
                } else {
                    throw new CDKException("Error while reading charge: expected double.");
                }
                IAtom atom = molecule.getAtom(atomCounter - 1);
                atom.setCharge(charge);
            }
        }
    }

    /**
     * Reads a set of vibrations into ChemModel.
     *
     * @param model the destination ChemModel
     * @throws java.io.IOException if an I/O error occurs
     */
//    private void readFrequencies(IChemModel model) throws IOException {
    /* This is yet to be ported. Vibrations don't exist yet in CDK.
   String line = input.readLine();
   line = input.readLine();
   line = input.readLine();
   line = input.readLine();
   line = input.readLine();
   while ((line != null) && line.startsWith(" Frequencies --")) {
       Vector currentVibs = new Vector();
       StringReader vibValRead = new StringReader(line.substring(15));
       StreamTokenizer token = new StreamTokenizer(vibValRead);
       while (token.nextToken() != StreamTokenizer.TT_EOF) {
           Vibration vib = new Vibration(Double.toString(token.nval));
           currentVibs.addElement(vib);
       }
       line = input.readLine(); // skip "Red. masses"
       line = input.readLine(); // skip "Rfc consts"
       line = input.readLine(); // skip "IR Inten"
       while (!line.startsWith(" Atom AN")) {
           // skip all lines upto and including the " Atom AN" line
           line = input.readLine(); // skip
       }
       for (int i = 0; i < frame.getAtomCount(); ++i) {
           line = input.readLine();
           StringReader vectorRead = new StringReader(line);
           token = new StreamTokenizer(vectorRead);
           token.nextToken();

           // ignore first token
           token.nextToken();

           // ignore second token
           for (int j = 0; j < currentVibs.size(); ++j) {
               double[] v = new double[3];
               if (token.nextToken() == StreamTokenizer.TT_NUMBER) {
                   v[0] = token.nval;
               } else {
                   throw new IOException("Error reading frequency");
               }
               if (token.nextToken() == StreamTokenizer.TT_NUMBER) {
                   v[1] = token.nval;
               } else {
                   throw new IOException("Error reading frequency");
               }
               if (token.nextToken() == StreamTokenizer.TT_NUMBER) {
                   v[2] = token.nval;
               } else {
                   throw new IOException("Error reading frequency");
               }
               ((Vibration) currentVibs.elementAt(j)).addAtomVector(v);
           }
       }
       for (int i = 0; i < currentVibs.size(); ++i) {
           frame.addVibration((Vibration) currentVibs.elementAt(i));
       }
       line = input.readLine();
       line = input.readLine();
       line = input.readLine();
   } */
//    }

    /**
     * Reads NMR nuclear shieldings.
     */
    private void readNMRData(IChemModel model, String labelLine) throws Exception {
        // FIXME: this is yet to be ported. CDK does not have shielding stuff.
        // Determine label for properties
        Pattern atomNumber = Pattern.compile("(\\d+)");
        Matcher matcher;
        int atomIndex = 0;
        IAtomContainer atomContainer = model.getMoleculeSet()
                .getAtomContainer(model.getMoleculeSet().getMoleculeCount() - 1);
        String line = input.readLine();
        boolean couplingTableFlag = false;
        List<Integer> couplingHeader=null;
        int matchEnd=0;
        // reset magnetic properties...
        // this is required for the l777 link
        // and it ensures that only the las set of parameters is recorded
        for (IAtom atom : atomContainer.atoms()){
            atom.setProperty(COUPLING_CONSTANTS, null);
            atom.setProperty(MAGNETIC_TENSOR, null);
        }
        // read the new magnetic parameters
        while (line != null && !REGEX_NMR_PARAMETERS_STOP.matcher(line).find()) {
            // extract the isotropic magnetic shielding tensor
            if (REGEX_MAGNETIC_TENSOR_TABLE.matcher(line).find()) {
                matcher = REGEX_MAGNETIC_TENSOR_TABLE.matcher(line);
                matcher.find();
                atomContainer.getAtom(Integer.parseInt(matcher.group(1)) - 1)
                        .setProperty(MAGNETIC_TENSOR, Double.parseDouble(
                                matcher.group(3).replace("D", "E")));
            }

            couplingTableFlag=(couplingTableFlag || REGEX_SS_COUPLING_FLAG.matcher(line).find());
            if(couplingTableFlag){
                // if it finds a line with 8 spaces, it must be the header of the coupling table
                // it is a bit weak but works for the moment
                if (REGEX_SS_COUPLING_HEADER_FLAG.matcher(line).find()){
                    couplingHeader = new ArrayList<Integer>();
                    matcher=atomNumber.matcher(line);
                    matchEnd=0;
                    // collect all the header values by matching a digit after the last position of the match
                    while (matcher.find(matchEnd)){
                        couplingHeader.add(Integer.parseInt(matcher.group(1))-1);
                        matchEnd=matcher.end();
                    }

                }
                // extract the coupling constants
                if(REGEX_SS_COUPLING_DIAGONAL.matcher(line).find()){
                    matcher = REGEX_SS_COUPLING_DIAGONAL.matcher(line);
                    matcher.find();
                    int atomID = Integer.parseInt(matcher.group(1))-1;
                    matchEnd=0;
                    int indexCount = 0;
                    // extract all the coupling values
                    while (REGEX_SS_COUPLING.matcher(line).find(matchEnd)){
                        matcher=REGEX_SS_COUPLING.matcher(line);
                        matcher.find(matchEnd);
                        matchEnd=matcher.end();
                        // skip when i==j
                        if(!couplingHeader.get(indexCount).equals(atomID)){
                            // define hashmap if property is not yet set
                            if(atomContainer.getAtom(atomID).getProperty(COUPLING_CONSTANTS)==null)
                                atomContainer.getAtom(atomID).setProperty(COUPLING_CONSTANTS,new HashMap<IAtom,Double>());
                            // add a mapping between atom and coupling value
                            ((HashMap<IAtom,Double>) atomContainer.getAtom(atomID).getProperty(COUPLING_CONSTANTS))
                                    .put(atomContainer.getAtom(couplingHeader.get(indexCount)),
                                            Double.parseDouble(matcher.group(1).replace("D", "E")));
                        }
                        indexCount++;
                    }
                }
            }
            line = input.readLine();
        }
    }

    /**
     * Select the theory and basis set from the first archive line.
     */
    private String parseLevelOfTheory(String line) {

        StringTokenizer st1 = new StringTokenizer(line, "\\");

        // Must contain at least 6 tokens
        if (st1.countTokens() < 6) {
            return null;
        }

        // Skip first four tokens
        for (int i = 0; i < 4; ++i) {
            st1.nextToken();
        }
        return st1.nextToken() + "/" + st1.nextToken();
    }


}
