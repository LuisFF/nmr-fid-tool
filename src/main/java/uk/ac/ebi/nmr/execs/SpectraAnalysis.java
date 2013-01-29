package uk.ac.ebi.nmr.execs;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.silent.ChemFile;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import uk.ac.ebi.nmr.io.G09OutputReader;
import uk.ac.ebi.nmr.io.L777OutputReader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 12/12/2012
 * Time: 16:38
 * To change this template use File | Settings | File Templates.
 */
public class SpectraAnalysis {

    public static void main(String[] args) throws Exception {


        if (args.length < 1)
            throw new Exception("Not enough arguments");
        if (!(new File(args[0]).isDirectory()))
            throw new Exception("The argument is not a directory");

        String workingDir = (new File(args[0])).getAbsolutePath();


        System.out.println("working directory: " + workingDir);

        FilenameFilter gaussianFilter = new FilenameFilter() {
            @Override
            public boolean accept(File file, String filename) {
                Pattern REGEX_GAUSSIAN_FILE = Pattern.compile("conf\\_\\d+(\\_solvent)?(-vibcorr)?\\.log");
                return REGEX_GAUSSIAN_FILE.matcher(filename).find();
            }
        };

        File[] gaussianFiles = (new File(workingDir)).listFiles(gaussianFilter);
        if (gaussianFiles.length == 0)
            throw new Exception("No gaussian file found");

        double[] lowestEnergy = new double[4];
        String[] lowestEnergyFile = new String[8];

        List<IAtomContainer> vacuumSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> solventSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> vacuumVibcorrSimulations = new ArrayList<IAtomContainer>();
        List<IAtomContainer> solventVibcorrSimulations = new ArrayList<IAtomContainer>();
        System.out.println(gaussianFiles.length);
        for (File file : gaussianFiles) {
            // check we are speaking about the gaussian files or not..
            if (!(file.getName().contains("err") || file.getName().contains("out"))) {
                BufferedReader bufferR = new BufferedReader(new FileReader(file));
                G09OutputReader g09OutputReader = new G09OutputReader(bufferR);
                L777OutputReader l777OutputReader = new L777OutputReader(bufferR);
                IChemFile chemFile = new ChemFile();
                // use a different buffer reader depeding on if this file an original guassian file or a l777 file
                if (file.getName().contains("vibcorr")) {
                    chemFile = l777OutputReader.read(chemFile);
                } else {
                    chemFile = g09OutputReader.read(chemFile);
                }
                // the reader will read a sequence of atomcontainer and only the last one will have all the information
                IChemModel model = ChemFileManipulator.getAllChemModels(chemFile)
                        .get(ChemFileManipulator.getAllChemModels(chemFile).size() - 1);

                if (file.getName().contains("vibcorr")) {
                    model.setProperty(G09OutputReader.NORMAL_TERMINATION,
                            (Boolean) model.getProperty(L777OutputReader.NORMAL_TERMINATION));
                    model.setProperty(G09OutputReader.STRUCTURE_ENERGY,
                            Double.parseDouble((String) model.getProperty(L777OutputReader.STRUCTURE_ENERGY)));

                }

                System.out.println("processing file: " + file.getName());
                System.out.println("Normal termination: " +
                        ((Boolean) model.getProperty(G09OutputReader.NORMAL_TERMINATION)));

                if ((Boolean) model.getProperty(G09OutputReader.NORMAL_TERMINATION)) {
//                    System.out.println("Energy of the geometry: " +
//                            model.getProperty(G09OutputReader.STRUCTURE_ENERGY));
                    //sort the atomcontainer
                    if (file.getName().contains("solvent")) {
                        if (file.getName().contains("vibcorr")) {
                            if (lowestEnergy[0] > (Double) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)) {
                                lowestEnergy[0] = (Double) model.getProperty(G09OutputReader.STRUCTURE_ENERGY);
                                lowestEnergyFile[0] = file.getName();
                            }
                        } else {
                            if (lowestEnergy[1] > Double.parseDouble((String) model
                                    .getProperty(G09OutputReader.STRUCTURE_ENERGY))) {
                                lowestEnergy[1] = Double.parseDouble((String) model
                                        .getProperty(G09OutputReader.STRUCTURE_ENERGY));
                                lowestEnergyFile[1] = file.getName();
                            }
                        }
                    } else {
                        if (file.getName().contains("vibcorr")) {
                            if (lowestEnergy[2] > (Double) model.getProperty(G09OutputReader.STRUCTURE_ENERGY)) {
                                lowestEnergy[2] = (Double) model.getProperty(G09OutputReader.STRUCTURE_ENERGY);
                                lowestEnergyFile[2] = file.getName();
                            }
                        } else {
                            if (lowestEnergy[3] > Double.parseDouble((String) model
                                    .getProperty(G09OutputReader.STRUCTURE_ENERGY))) {
                                lowestEnergy[3] = Double.parseDouble((String) model
                                        .getProperty(G09OutputReader.STRUCTURE_ENERGY));
                                lowestEnergyFile[3] = file.getName();
                            }
                        }
                    }
                }
            } // write down the NMR data
//            if(!new File(workingDir+"/"+ NMRDBSimulatorAllConf.SimulationFiles.values()).exists()) {
//                System.err.println("Spectra file \" + filename + \" does not exist");
//            }

        }
        int count=4;
        for (NMRDBSimulatorAllConf.SimulationFiles files: NMRDBSimulatorAllConf.SimulationFiles.values()){
            if (!new File(workingDir + "/" + files.filename()).exists()) {
                System.err.println("Spectra file " + files.filename()+ " does not exist");
            } else {
                lowestEnergyFile[count]=files.filename();
                count++;
            }
        }

        writeMatLabFile(lowestEnergyFile, workingDir);


//        if(vacuumSimulations.size()==0)
//            throw new Exception("No simulations in vacuum conditions found");
//        if(solventSimulations.size()==0)
//            throw new Exception("No simulations in solvent conditions found");
//        if(solventVibcorrSimulations.size()==0)
//            throw new Exception("No simulations in solvent and vibrational corrections conditions found");
//        if(vacuumVibcorrSimulations.size()==0)
//            throw new Exception("No simulations in vacuum and vibrational corrections conditions found");


    }

    private static void writeMatLabFile(String[] lowestEnergyFile, String workingDir) {
        FileWriter fstream = null;
        try {
            fstream = new FileWriter(workingDir + "/analysisWCC.m");
            System.out.println("bu");

            BufferedWriter out = new BufferedWriter(fstream);
            int count=1;
            out.write("expSpectrum=csvread('experimental-spectrum.csv');\n");
            out.write("shifts=expSpectrum((expSpectrum(:,1)<9.8 & expSpectrum(:,1)>0.2),1);\n"+
                      "experimental=normalizeIntensity(expSpectrum((expSpectrum(:,1)<9.8 & expSpectrum(:,1)>0.2),3));\n");
            String spectraLabelMap="spectraLable={};\n";
            String normalizeSpectra="";
            String wccAnalysis="";
//            String figure="figure;\nplot(shifts,experimental";
            String figure="";
//            String legend="legend('experimental'";
//            String figure="figure;\nsubplot(XX, 1, "+count+");\nplot(shifts,experimental);\n"+
//                    "set(gca,'XDir','Reverse');\nxlabel('ppm');\nset(gca,'Ticklength',[0 0]);\n"+
//                    "axis([0 9 -0.05 1.2]);\nset(gca, 'LineWidth', 1);\ntitle('Experimental');\n";
            for (String filename : lowestEnergyFile) {
                filename = filename.replaceFirst("-sim.txt", "-spectra.txt");
                filename = filename.replaceFirst(".log", "-spectra.txt");
                String spectraName = filename.replace("-spectra.txt","");
                spectraName=spectraName.replaceAll("-","_");
                System.out.println(filename);
                if (!new File(workingDir + "/" + filename).exists()) {
                    System.err.println("Spectra file " + filename + " does not exist");
                } else {
                    out.write(spectraName+"=csvread('"+ filename + "');\n");
                    spectraLabelMap+="spectraLable{"+count+"}='"+spectraName+"';\n";
                    normalizeSpectra+="simulated(:,"+count+")=normalizeIntensity("+
                            "interp1("+spectraName+"(:,1),"+spectraName+"(:,2),shifts,'cubic'));\n";
                    wccAnalysis+="disp(sprintf('WCC for "+spectraName+
                            ": %d\n', wccAnalysis02(1,0.01,shifts', experimental,simulated(:,"+count+ ")')));\n";
//                    figure+="subplot(XX, 1, "+(count+1)+");\nplot(shifts,simulated(:,"+count+"));\n"+
                    figure+="figure;\nset(gcf,'numbertitle','off','name','"+spectraName+"');\n"+
                            "plot(shifts,experimental,shifts,simulated(:,"+count+"));\n"+
                            "set(gca,'XDir','Reverse');\nxlabel('ppm');\nset(gca,'Ticklength',[0 0]);\n"+
                            "axis([0 9 -0.05 1.2]);\nset(gca, 'LineWidth', 1);\ntitle(regexprep('"+ spectraName+
                            "', '_', '\\\\_'));\n" +
                            "fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 8.5]);\n"+
                            "print(gcf, '-dpdf', '-r300', '"+spectraName+".pdf');\n";
//                    figure+=",shifts,simulated(:,"+count+")";
//                    legend+=","+spectraName;
                    count++;
                }
            }
//            figure=figure.replaceAll("XX",String.valueOf(count));
//            figure+=");\nset(gca,'XDir','Reverse');\nxlabel('ppm');\nset(gca,'Ticklength',[0 0]);\n"+
//                    "axis([0 9 -0.05 1.2]);\nset(gca, 'LineWidth', 1);\n";
//            legend+=",'Location','NorthEast');\n";
            out.write(spectraLabelMap+normalizeSpectra+figure);
            // plot just the regions around the peaks
            out.write("selection=experimental>0.05; \n" +
                    "for n=1:size(simulated,2)\n" +
                    "selection2=simulated(:,n)>0.05;\n" +
                    "selection2=selection| selection2';\n" +
                    "count=1;\n" +
                    "status=false;\n" +
                    "finalSelection(1:length(selection2))=false;\n" +
                    "for i=1:length(selection2)     \n" +
                    "if((xor(status,selection2(i))))\n" +
                    "status=selection2(i);\n" +
                    "if(selection2(i))\n" +
                    "limits(count,1)=max(1,i-250);\n" +
                    "else\n" +
                    "limits(count,2)=min(length(selection2),i+250);\n" +
                    "count=count+1;\n" +
                    "end\n" +
                    "end\n" +
                    "end\n" +
                    "\n" +
                    "% remove overlap between intervals\n" +
                    "limitsupdate=limits(1,:);\n" +
                    "for i=2:length(limits)\n" +
                    "if(limits(i,1)<limitsupdate(size(limitsupdate,1),2))\n" +
                    "limitsupdate(size(limitsupdate,1),2)=limits(i,2);\n" +
                    "%limits(1:i-1,:)=limits(1:i-1,:);\n" +
                    "%limits(i:(length(limits)-1),:)=limits((i+1):length(limits),:);\n" +
                    "%limits=limits(1:(length(limits)-1),:);\n" +
                    "else\n" +
                    "limitsupdate(size(limitsupdate,1)+1,:)=limits(i,:);\n" +
                    "end\n" +
                    "end\n" +
                    "\n" +
                    "%subplot....\n" +
                    "\n" +
                    "figure;\n" +
                    "set(gcf,'numbertitle','off','name',char(spectraLable(n))) \n" +
                    "for i=1:size(limitsupdate,1)\n" +
                    "subplot(1,size(limitsupdate,1), i);\n" +
                    "plot(shifts(limitsupdate(i,1):limitsupdate(i,2)),experimental(limitsupdate(i,1):limitsupdate(i,2)),"+
                    "shifts(limitsupdate(i,1):limitsupdate(i,2)),simulated((limitsupdate(i,1):limitsupdate(i,2)),n));\n" +
                    "set(gca,'XDir','Reverse');\n" +
                    "xlabel('ppm');\n" +
                    "set(gca,'Ticklength',[0 0]);\n" +
                    "set(gca, 'LineWidth', 1);\n" +
                    "ylim([-0.05 1.2]);\n"+
                    "fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 8.5]);\n"+
                    "print(gcf, '-dpdf', '-r300', strcat(char(spectraLable(n)),'-zoom','.pdf'));\n" +
                    "end\n" +
                    "end\n"
            );
            out.close();

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
