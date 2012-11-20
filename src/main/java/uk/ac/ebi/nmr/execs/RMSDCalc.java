/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.ebi.nmr.execs;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.geometry.alignment.KabschAlignment;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import uk.ac.ebi.nmr.io.G09InputReader;
import uk.ac.ebi.nmr.io.LasCMDReader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author ldpf
 */
public class RMSDCalc {


    public static void main(String[] args) throws Exception {
        FileInputStream iStream = null;
        Pattern energyPattern = Pattern.compile("Energy\\(kCal\\/mol\\) = (-?\\d+\\.\\d+)");
        Matcher matcher = null;
        Integer[] lowestEnergyIndex2 = null;
        if (args.length < 2)
            throw new Exception("Not enough arguments");

        int kclusters = Integer.parseInt(args[1]);
        String filepath = args[0];
        // one has to work with the full path otherwise it is not possible to extract the working directory
        String workingDir = (new File((new File(filepath)).getAbsolutePath())).getParent();


        try {

            System.out.println("Reading file: " + args[0]);
            System.out.println("working directory" + (new File(workingDir)).getAbsolutePath());
            System.out.println("Nb of clusters: " + args[1]);
            IChemFile chemFile = new ChemFile();
            LasCMDReader lasCMDReader = new LasCMDReader(new FileReader(filepath));
            chemFile = (IChemFile) lasCMDReader.read(chemFile);
            IAtomContainer atomContainerOriginal = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                    ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1);
            lasCMDReader.close();

            String inchiOriginal = InChIGeneratorFactory.getInstance()
                    .getInChIGenerator(atomContainerOriginal).getInchi();

            FilenameFilter gaussianFilter = new FilenameFilter() {
                @Override
                public boolean accept(File file, String filename) {

                    return filename.endsWith("com") && filename.contains("geoms");
                }
            };

            File[] gaussianFiles = (new File(workingDir)).listFiles(gaussianFilter);
            if(gaussianFiles.length==0)
                throw new Exception("No geometry file found");
            System.out.println("Number of gaussian files: "+gaussianFiles.length );
//            iStream = new FileInputStream(new File(filepath));
//            MDLReader mdlReader = new MDLReader(iStream);
            G09InputReader g09InputReader = null;
            List<IAtomContainer> atomContainers = new ArrayList<IAtomContainer>(gaussianFiles.length);
            List<Double> energies = new ArrayList<Double>(gaussianFiles.length);
            /*
            extract information from the gaussian files
             */
            for (File file : gaussianFiles) {
                g09InputReader = new G09InputReader(new FileReader(file));
                chemFile = (IChemFile) g09InputReader.read(chemFile);
                IAtomContainer conformer = ChemFileManipulator.getAllAtomContainers(chemFile).get(
                        ChemFileManipulator.getAllAtomContainers(chemFile).size() - 1
                );
                g09InputReader.close();
                String header = (String) ChemFileManipulator.getAllChemModels(chemFile)
                        .get(ChemFileManipulator.getAllChemModels(chemFile).size() - 1)
                        .getProperty(G09InputReader.HEADER);

                matcher = energyPattern.matcher(header);

                if (!matcher.find())
                    throw new Exception("No energy found in the header of Gaussian");
                if (!InChIGeneratorFactory.getInstance()
                        .getInChIGenerator(conformer).getInchi().equals(inchiOriginal))
                    System.err.println("[WARN] Conformer structure changed from the original");
//                System.out.println(smiler.createSMILES(conformer));
                atomContainers.add(conformer);
                energies.add(Double.parseDouble(matcher.group(1)));


            }

            System.out.println("Number of geometries: "+atomContainers.size());
            System.out.println("InChI for the first geomtry: " + InChIGeneratorFactory.getInstance()
                    .getInChIGenerator(atomContainers.get(0)).getInchi());
            /*
            Skip the current implementation of the clustering because it's not working for large data :S
             */
//            RMSDCluster rmsdCluster = new RMSDCluster(atomContainers);
//            rmsdCluster.clusterRMSDs(kclusters);

            /*
            Using the implementation of the method by Park and Jun
            It is required to calculate the RMSDs
             */
//            Integer[] clusterAsignment = rmsdCluster.getClusterAssignment();
            // perform the clustering only if one has less than 7 geometries
            if (atomContainers.size() > 7) {


                // get the molecule with the lowest energy for each cluster
                Double[] lowestE = new Double[kclusters];
                IAtomContainer[] lowestEAtomContainer = new IAtomContainer[kclusters];
                Integer[] lowestEnergyIndex = new Integer[kclusters];

//            for (int i = 0; i < atomContainers.size(); i++) {
//                if (lowestE[clusterAsignment[i] - 1] == null)
//                    lowestE[clusterAsignment[i] - 1] = 1000.0;
//                // I had to store the energy as the name because Obabel does not allow to collect this info
//                if (lowestE[clusterAsignment[i] - 1] > Double.parseDouble((String) atomContainers.get(i).getProperty("cdk:Title"))) {
//                    lowestE[clusterAsignment[i] - 1] = Double.parseDouble((String) atomContainers.get(i).getProperty("cdk:Title"));
//                    lowestEAtomContainer[clusterAsignment[i] - 1] = atomContainers.get(i);
//                    lowestEnergyIndex[clusterAsignment[i] - 1] = i;
//                }
//            }
//
////            for (int i = 0 ; i < kclusters; i++){
////                // I had to store the energy as the name because Obabel does not allow to collect this info
////                System.out.println("cluster "+(i+1)+" : "+lowestE[i]);
////            }
//            for (int i = 0; i < kclusters; i++) {
//                // I had to store the energy as the name because Obabel does not allow to collect this info
//                System.out.println("cluster " + (i + 1) + " : " + lowestE[i]);
//            }
//
//            System.out.print("medoids: " +
//                    rmsdCluster.getMedoidIndexes()[0] +
//                    "   " +
//                    rmsdCluster.getMedoidIndexes()[1] +
//                    "   " +
//                    rmsdCluster.getMedoidIndexes()[2] +
//                    "   " +
//                    rmsdCluster.getMedoidIndexes()[3] +
//                    "   " +
//                    rmsdCluster.getMedoidIndexes()[4] + "\n");

                /*
               do the clustering using R and the RMSDs obtained from the RMSDCluster
                */
                // calculate the RMSDs
                double[][] rmsdMatrix = calculateRMSDs(atomContainers);
                // Create file dissimilarity file
                writeDisMatrix(workingDir, rmsdMatrix, energies);
                // run PAM in R
                runRClustering(workingDir, kclusters);
                //parse the R output
                Integer[] assignment = parseRClustering(atomContainers.size(), workingDir);

                // obtain the lowest energy atomcontainers
                Double[] lowestEnergy2 = new Double[kclusters];

                lowestEnergyIndex2 = new Integer[kclusters];
                for (int i = 0; i < assignment.length; i++) {
                    if (lowestEnergy2[assignment[i] - 1] == null)
                        lowestEnergy2[assignment[i] - 1] = 1000.0;
                    if (lowestEnergy2[assignment[i] - 1] > energies.get(i)) {
                        lowestEnergy2[assignment[i] - 1] = energies.get(i);
                        lowestEnergyIndex2[assignment[i] - 1] = i;
                    }

                }
//            for (int i = 0; i < kclusters; i++) {
//                System.out.println(lowestEnergyIndex[i] + " " + lowestE[i] + " | " + lowestEnergyIndex2[i] + " " + lowestEnergy2[i]);
//                System.out.println(rmsdCluster.getMedoidIndexes()[i] + " " + assignment[rmsdCluster.getMedoidIndexes()[i]]);
//            }
            } else {
                lowestEnergyIndex2 = new Integer[atomContainers.size()];
                for (int i = 0; i < lowestEnergyIndex2.length; i++) {
                    lowestEnergyIndex2[i] = i;
                }
            }

            writeGaussian(workingDir, atomContainers, lowestEnergyIndex2, energies);

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (CDKException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }

    private static double[][] calculateRMSDs(List<IAtomContainer> atomContainers) throws Exception {
        double[][] rmsds = new double[atomContainers.size()][atomContainers.size()];

        for (int i = 0; i < atomContainers.size() - 1; i++) {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainers.get(i));
            CDKHueckelAromaticityDetector.detectAromaticity(atomContainers.get(i));
            if (!GeometryTools.has3DCoordinates(atomContainers.get(i))) {
                throw new Exception("Molecule missing 3D coordinates");
            }
            for (int j = i + 1; j < atomContainers.size(); j++) {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainers.get(j));
                CDKHueckelAromaticityDetector.detectAromaticity(atomContainers.get(j));
                if (!GeometryTools.has3DCoordinates(atomContainers.get(j))) {
                    throw new Exception("Molecule missing 3D coordinates");
                }
                /*
                Isomorphic alignment
                 */
//                    Map<Integer, Integer> map = new HashMap<Integer, Integer>();
//                    if (atomContainers.get(i).getAtomCount() > atomContainers.get(j).getAtomCount()) {
//                        map = AtomMappingTools
//                                .mapAtomsOfAlignedStructures(
//                                        atomContainers.get(i), atomContainers.get(j), map);
//
//                        rmsds[i][j] = GeometryTools.getAllAtomRMSD(
//                                atomContainers.get(i),
//                                atomContainers.get(j),
//                                map,
//                                true);
//                        rmsds[j][i] = rmsds[i][j];
//
//                    } else {
//                        map = AtomMappingTools
//                                .mapAtomsOfAlignedStructures(
//                                        atomContainers.get(j), atomContainers.get(i), map);
//                        rmsds[i][j] = GeometryTools.getAllAtomRMSD(
//                                atomContainers.get(j),
//                                atomContainers.get(i),
//                                map,
//                                true);
//                        rmsds[j][i]=rmsds[i][j];
//                    }

                /*
                for KabschAlignment
                 */
                KabschAlignment sa = new KabschAlignment(atomContainers.get(i), atomContainers.get(j));
                sa.align();
                rmsds[i][j] = sa.getRMSD();
                rmsds[j][i] = rmsds[i][j];
//                    if(rmsds[i][j]){
//                        throw new Exception("Problem aligning molecules");
//                    }
            }
        }
        return rmsds;  //To change body of created methods use File | Settings | File Templates.
    }

    private static void writeGaussian(String workingDir,
                                      List<IAtomContainer> atomContainers,
                                      Integer[] lowestEnergyIndex2,
                                      List<Double> energies) throws IOException, CDKException {

        for (int index : lowestEnergyIndex2) {
            //write the molfile
            MDLV2000Writer molWriter = new MDLV2000Writer(new FileWriter(
                    new File(workingDir + "/conf_" + index + ".mol")));
            molWriter.write(atomContainers.get(index));
            molWriter.close();

            //
            String header1 = "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + ".chk\n" +
                    "# B3LYP/6-31+G Opt=loose int=grid=sg1 Test\n\n" +
                    "Geometry generated by LasCMD Energy (kcal/mol) = " +
                    energies.get(index)+
                    "\n[Gas-Phase]\n" +
                    "\n" +
                    "0,1\n";
            String header2 = "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + "_solvent.chk\n" +
                    "# B3LYP/6-31+G Opt=loose int=grid=sg1 Test SCRF=(PCM,solvent=Water)\n\n" +
                    "Geometry generated by LasCMD Energy (kcal/mol) = " +
                    energies.get(index) +
                    "\n[Water]\n" +
                    "\n" +
                    "0,1\n";

            String footer1 = "\n\n--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%nosave\n" +
                    "%chk=conf_" + index + ".chk\n" +
                    "# B3LYP/6-311++G** Opt=tight int=ultrafinegrid Test geom=allcheck guess=read\n" +
                    "\n" +
                    "\n" +
                    " perform a better optimization\n" +
                    "\n" +
                    "--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + ".chk\n" +
                    "# B3LYP/6-311++G** Freq Test Geom=allcheck Guess=Tcheck\n" +
                    "\n" +
                    " calculate Frequencies and Energies\n" +
                    "\n" +
                    "--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + ".chk\n" +
                    "# B3LYP/6-311++G** NMR Test Geom=allcheck Guess=read\n" +
                    "\n" +
                    " calculate the Magnetic tensors\n" +
                    "\n" +
                    "--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + ".chk\n" +
                    "# B3LYP/6-311++G** NMR=mixed Test Geom=allcheck Guess=read\n" +
                    "\n" +
                    " recalculate the Magnetic tensors and the spin-spin couplings\n" +
                    "\n" +
                    "\n";

            String footer2 = "\n\n--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%nosave\n" +
                    "%chk=conf_" + index + "_solvent.chk\n" +
                    "# B3LYP/6-311++G** Opt=tight int=ultrafinegrid Test geom=allcheck guess=read SCRF=(PCM,solvent=Water)\n" +
                    "\n" +
                    "\n" +
                    " perform a better optimization\n" +
                    "\n" +
                    "--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + "_solvent.chk\n" +
                    "# B3LYP/6-311++G** Freq Test Geom=allcheck Guess=Tcheck SCRF=(PCM,solvent=Water)\n" +
                    "\n" +
                    " calculate Frequencies and Energies\n" +
                    "\n" +
                    "--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + "_solvent.chk\n" +
                    "# B3LYP/6-311++G** NMR Test Geom=allcheck Guess=read SCRF=(PCM,solvent=Water)\n" +
                    "\n" +
                    " calculate the Magnetic tensors\n" +
                    "\n" +
                    "--link1--\n" +
                    "%mem=64GB\n" +
                    "%nproc=12\n" +
                    "%chk=conf_" + index + "_solvent.chk\n" +
                    "# B3LYP/6-311++G** NMR=mixed Test Geom=allcheck Guess=read SCRF=(PCM,solvent=Water)\n" +
                    "\n" +
                    " recalculate the Magnetic tensors and the spin-spin couplings\n" +
                    "\n" +
                    "\n";

            FileWriter fstream1 = new FileWriter(workingDir + "/conf_" + index + ".inp");
            FileWriter fstream2 = new FileWriter(workingDir + "/conf_" + index + "_solvent.inp");
            BufferedWriter out1 = new BufferedWriter(fstream1);
            BufferedWriter out2 = new BufferedWriter(fstream2);
            out1.write(header1);
            out2.write(header2);
            String atomCoordinates = "";
            for (int i = 0; i < atomContainers.get(index).getAtomCount(); i++) {
                atomCoordinates += atomContainers.get(index).getAtom(i).getSymbol() + "," +
                        atomContainers.get(index).getAtom(i).getPoint3d().x + "," +
                        atomContainers.get(index).getAtom(i).getPoint3d().y + "," +
                        atomContainers.get(index).getAtom(i).getPoint3d().z + "\n";
            }

            out1.write(atomCoordinates);
            out2.write(atomCoordinates);
            out1.write(footer1);
            out2.write(footer2);
            out1.close();
            out2.close();
        }
    }

    private static void writeDisMatrix(String workingDir, double[][] rmsds, List<Double> energies) throws IOException {
        FileWriter fstream = new FileWriter(workingDir + "/dissimilarity.txt");
        FileWriter fstreamEnergies = new FileWriter(workingDir + "/energies.txt");

        BufferedWriter out = new BufferedWriter(fstream);
        BufferedWriter out2 = new BufferedWriter(fstreamEnergies);
        PrintWriter printWriter = new PrintWriter(out);
        PrintWriter printWriter2 = new PrintWriter(out2);
        for (int i = 0; i < rmsds.length; i++) {
            printWriter2.printf("%.6f%n", energies.get(i));

            for (int j = 0; j < rmsds.length; j++) {
//                if (rmsds[i][j] == null)
//                      printWriter.printf("%.7f",0);
//                    out.write("0\t");
                if (j > rmsds.length - 1)
                    printWriter.printf("%.7f", rmsds[i][j]);
//                    out.write(rmsdCluster.getRmsds()[i][j].toString());
                else
                    printWriter.printf("%.7f\t", rmsds[i][j]);
//                    out.write(rmsdCluster.getRmsds()[i][j] + "\t");
            }
            printWriter.printf("%n", 0);
//            out.write("\n");
        }
        out2.close();
        out.close();
    }

    /**
     * parse the results of the clustering and assign each conformation the corresponding cluster.
     *
     * @param size
     * @param workingDir
     * @return
     * @throws java.io.IOException
     */
    private static Integer[] parseRClustering(int size, String workingDir) throws IOException {

        FileReader reader = new FileReader(workingDir + "/clusterResults.txt");
        BufferedReader buffReader = new BufferedReader(reader);
        String line = "";
        Integer[] assignment = new Integer[size];
        Pattern p = Pattern.compile("^\"V([0-9]+)\"");
        while ((line = buffReader.readLine()) != null) {

            if (line.startsWith("\"V")) {
                String[] vector = line.split(" ");
                Matcher m = p.matcher(vector[0]);
                if (m.find()) {
                    assignment[Integer.parseInt(m.group(1)) - 1] = Integer.parseInt(vector[1]);
                }
            }
        }
        buffReader.close();
        return assignment;
    }

    private static void runRClustering(String workingDir, int clusterNb) throws IOException, InterruptedException {
        FileWriter fstream = new FileWriter(workingDir + "/pam_clustering.R");
        BufferedWriter out = new BufferedWriter(fstream);
        out.write("library(pamr)\n" +
                "library(StatDA)\n" +
                "data<-read.table(\"" + workingDir + "/dissimilarity.txt\")\n" +
                "energies<-read.table(\"" + workingDir + "/energies.txt\")\n" +
                "clusters <- pam(as.dist(data), " + clusterNb + ", diss=1)\n" +
                "result <- clusters[3]\n" +
                "postscript(\"" + workingDir + "/pam_clustering.eps\")\n" +
                "clusplot(clusters)\n" +
                "dev.off()\n" +
                "postscript(\"" + workingDir + "/pam_clustering_energies_1.eps\")\n" +
                "edaplot(energies[stack(clusters[3])[1]==1])\n" +
                "dev.off()\n" +
                "postscript(\"" + workingDir + "/pam_clustering_energies_2.eps\")\n" +
                "edaplot(energies[stack(clusters[3])[1]==2])\n" +
                "dev.off()\n" +
                "postscript(\"" + workingDir + "/pam_clustering_energies_3.eps\")\n" +
                "edaplot(energies[stack(clusters[3])[1]==3])\n" +
                "dev.off()\n" +
                "postscript(\"" + workingDir + "/pam_clustering_energies_4.eps\")\n" +
                "edaplot(energies[stack(clusters[3])[1]==4])\n" +
                "dev.off()\n" +
                "postscript(\"" + workingDir + "/pam_clustering_energies_5.eps\")\n" +
                "edaplot(energies[stack(clusters[3])[1]==5])\n" +
                "dev.off()\n" +
                "postscript(\"" + workingDir + "/energy_distribution.eps\")\n" +
                "edaplot(energies[,1])\n" +
                "dev.off()\n" +
                "\n" +
                "#Sys.sleep(10)\n" +
                "write.table(result, file = \"" + workingDir + "/clusterResults.txt\")\n" +
                "q(save = \"yes\", status = 0, runLast = TRUE)");
        out.close();
        // run R to get the clustering
        Process p = Runtime.getRuntime().exec("R CMD BATCH " + workingDir + "/pam_clustering.R");
        //wai to finish the R
        p.waitFor();
    }

}
