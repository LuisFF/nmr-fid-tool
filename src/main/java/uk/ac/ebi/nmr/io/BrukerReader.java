package uk.ac.ebi.nmr.io;



import java.io.File;
import java.io.FileNotFoundException;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 21/09/2012
 * Time: 14:46
 * To change this template use File | Settings | File Templates.
 */
public class BrukerReader {

    private File fileFID = null;
    private File fileACQU = null;
    private File filePROC = null;

    //TODO review regexp

    // parameters from acqu
    private final static Pattern REGEXP_SFO1 = Pattern.compile("\\#\\#\\$SFO1= (-?\\d+\\.\\d+)"); //transmitter frequency
    private final static Pattern REGEXP_SFO2 = Pattern.compile("\\#\\#\\$SFO2= (-?\\d+\\.\\d+)"); //decoupler frequency
    private final static Pattern REGEXP_SFO3 = Pattern.compile("\\#\\#\\$SFO3= (\\d+\\.\\d+)"); //second decoupler frequency
    private final static Pattern REGEXP_O1 = Pattern.compile("\\#\\#\\$O1= (\\d+\\.\\d+)"); //frequency offset in Hz
    private final static Pattern REGEXP_SW = Pattern.compile("\\#\\#\\$SW= (\\d+\\.\\d+)"); //spectral width (ppm)
    private final static Pattern REGEXP_TD = Pattern.compile("\\#\\#\\$TD= (\\d+)"); //acquired points (real+imaginary)
    private final static Pattern REGEXP_DECIM = Pattern.compile("\\#\\#\\$DECIM= (-?\\d+)"); //DSP decimation factor
    private final static Pattern REGEXP_DSPFVS = Pattern.compile("\\#\\#\\$DSPFVS= (-?\\d+)"); //DSP firmware version
    private final static Pattern REGEXP_GRPDLY = Pattern.compile("\\#\\#\\$GRPDLY= (-?\\d+)"); //DSP group delay
    private final static Pattern REGEXP_BYTORDA = Pattern.compile("\\#\\#\\$BYTORDA= (\\d+)"); //byte order
    // variables not yet defined in Experiment
    private final static Pattern REGEXP_AQ_MODE = Pattern.compile("\\#\\#\\$AQ\\_mod= (\\d+)"); //acquisition mode
    private final static Pattern REGEXP_DIGMOD = Pattern.compile("\\#\\#\\$DIGMOD= (\\d+)"); //filter type
    private final static Pattern REGEXP_NS = Pattern.compile("\\#\\#\\$NS= (\\d+)"); //number of scans
    //TODO review REGEXP_PULPROG
    // examples of REGEXP_PULPROG : <zg> <cosydfph> <bs_hsqcetgpsi>; basically a word between < >
    private final static Pattern REGEXP_PULPROG = Pattern.compile("\\#\\#\\$PULPROG= "); //pulse program
    //TODO review REGEXP_NUC1
    // examples of REGEXP_NUC1 : <1H>; basically <isotope number + element>
    private final static Pattern REGEXP_NUC1 = Pattern.compile("\\#\\#\\$NUC1= "); // observed nucleus
    //TODO review REGEXP_INSTRUM
    // examples of REGEXP_INSTRUM : <amx500> ; basically <machine name>
    private final static Pattern REGEXP_INSTRUM = Pattern.compile("\\#\\#\\$INSTRUM= "); // instrument name
    private final static Pattern REGEXP_DTYPA = Pattern.compile("\\#\\#\\$DTYPA= (\\d+)"); //data type (0 -> 32 bit int, 1 -> 64 bit double)
    //TODO review REGEXP_SOLVENT
    // examples of REGEXP_SOLVENT : <DMSO> ; basically <solvent name>
    private final static Pattern REGEXP_SOLVENT = Pattern.compile("\\#\\#\\$SOLVENT= "); // instrument name
    //TODO review REGEXP_PROBHD
    // examples of REGEXP_PROBHD : <32> <>; basically <digit?>
    private final static Pattern REGEXP_PROBHD = Pattern.compile("\\#\\#\\$PROBHD= "); // probehead
    //TODO review REGEXP_ORIGIN
    // examples of REGEXP_ORIGIN : Bruker Analytik GmbH; basically a name
    private final static Pattern REGEXP_ORIGIN = Pattern.compile("\\#\\#\\$ORIGIN= "); // origin
    //TODO review REGEXP_OWNER
    // examples of REGEXP_OWNER : guest; basically the used ID
    private final static Pattern REGEXP_OWNER = Pattern.compile("\\#\\#\\$OWNER= "); // owner



    // parameters from proc
    private final static Pattern REGEXP_SI = Pattern.compile("\\#\\#\\$SI= (\\d+\\.\\d+)"); //transform size (complex)
    private final static Pattern REGEXP_SF = Pattern.compile("\\#\\#\\$SF= (\\d+\\.\\d+)"); //frequency of 0 ppm
    private final static Pattern REGEXP_GB = Pattern.compile("\\#\\#\\$GB= (\\d+\\.\\d+)"); //GB-factor (Gain?)
    private final static Pattern REGEXP_LB = Pattern.compile("\\#\\#\\$LB= (\\d+)"); //line broadening
    // variables not yet defined in Experiment
    private final static Pattern REGEXP_WDW = Pattern.compile("\\#\\#\\$WDW= (\\d+)"); //window function type
    private final static Pattern REGEXP_PH_MODE = Pattern.compile("\\#\\#\\$PH\\_mod= (\\d+)"); //phasing type
    private final static Pattern REGEXP_PHC0 = Pattern.compile("\\#\\#\\$PHC0= (-?\\d+\\.\\d+)"); //zero order phase
    private final static Pattern REGEXP_PHC1 = Pattern.compile("\\#\\#\\$PHC1= (-?\\d+\\.\\d+)"); //first order phase
    private final static Pattern REGEXP_SSB = Pattern.compile("\\#\\#\\$SSB= (-?\\d+\\.\\d+)"); //sine bell shift
    private final static Pattern REGEXP_MC2 = Pattern.compile("\\#\\#\\$MC2= (\\d+)"); //F1 detection mode



    private Experiment experiment;
    private Process process;





    public BrukerReader() {
    }

    public BrukerReader(String filename) throws FileNotFoundException {

        this.fileFID = new File(filename);
        String workingDIR = fileFID.getParent();
        this.fileACQU = new File(workingDIR+"acqu");
        // TODO make other PROC files available
        this.filePROC = new File(workingDIR+"/pdata/1/proc");

    }

    public void read(){
        this.experiment = new Experiment();
        this.process = new Process();

        readACQU();
        // for the moment set only one processing though it is irrelevant for now
        readPROC(1);
    }

    /**
     * read the acq file and extract the parameters
     */
    private void readACQU () {

    }

    /**
     * read proc file and extract the paramenters. Note that can be more than one processing.
     * @param processingNb
     */
    private void readPROC (int processingNb) {

    }



    private class Process {
        private int transformSize;
        private double originFreq;
        private double gbFactor;
        private double lineBroadening;

        private Process() {
        }
    }

    private class Experiment {

        private double transmitterFreq;
        private double decouplerFreq;
        private double secondDecouplerFreq;
        private double offsetFreq;
        private double spectralWidth;



        private int acquiredPoints;
        private int dspDecimationFactor;
        private int dspFirmwareVersion;
        private double dspGroupDelay;

        private int byteOrder;


        private Experiment() {
        }
    }

}
