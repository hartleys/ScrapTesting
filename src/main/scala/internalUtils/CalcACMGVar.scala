package internalUtils

import java.io.BufferedReader;
import java.io.BufferedInputStream;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.Writer;
import java.io.File;
import java.util.zip.GZIPOutputStream;
import java.util.zip.GZIPInputStream;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;

import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import scala.collection.JavaConversions._
import internalUtils.optionHolder._;
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;

import htsjdk.variant._;
import htsjdk.variant.variantcontext._;
import htsjdk.variant.vcf._;

import VcfTool._;

import internalUtils.commandLineUI._;

object CalcACMGVar {
  
  class CmdAssessACMG extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "CmdAssessACMG", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "BETA: This function consolidates information from a wide variety of different input files and attempts to calculate a "+
                        "subset of the ACMG guidelines criteria. It then attempts to assign pathogenicity scores. " + BETA_WARNING,   
          argList = 
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "clinVarVcf", 
                                         arg = List("--clinVarVcf"), 
                                         valueName = "clinVarVcf.vcf",  
                                         argDesc =  "Processed clinvar variant vcf file. This file must have been processed by the addTxInfoToVCF command."
                                        ) ::  
                    new BinaryOptionArgument[String](
                                         name = "txToGeneFile", 
                                         arg = List("--txToGeneFile"), 
                                         valueName = "txToGene.txt",  
                                         argDesc =  "File containing the mapping of transcript names to gene symbols. This file must have 2 columns: the txID and the geneID. No header line."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "rmskFile", 
                                         arg = List("--rmskFile"), 
                                         valueName = "rmsk.txt.gz",  
                                         argDesc =  "rmsk.txt.gz file, from UCSC."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "toleranceFile", 
                                         arg = List("--toleranceFile"), 
                                         valueName = "toleranceFile.txt",  
                                         argDesc =  "This file must contain three columns (labelled in a header line): geneID (gene symbol), LOFtolerant, and MIStolerant. Genes that are not included in this list will be assumed to be non-tolerant."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "domainFile", 
                                         arg = List("--domainFile"), 
                                         valueName = "domainFile.txt",  
                                         argDesc =  "This file must contain at least four columns (labelled in a header line): chrom, start, end, and domainID."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "conservedElementFile", 
                                         arg = List("--conservedElementFile"), 
                                         valueName = "conservedElementFile.txt",  
                                         argDesc =  "This file contains the spans for the conserved element regions found by GERP. This file must contain 3 columns (no header line): chrom, start end."
                                        ) :: 
                    new BinaryArgument[List[String]](name = "ctrlAlleFreqKeys",
                                           arg = List("--ctrlAlleFreqKeys"),  
                                           valueName = "key1,key2,...", 
                                           argDesc = "List of VCF INFO tags tcontaining the allele frequencies for the control datasets.",
                                           defaultValue = Some(List("1KG_AF","ESP_EA_AF","ExAC_ALL"))
                                           ) ::
                    new BinaryArgument[Double](name = "BA1_AF",
                                           arg = List("--BA1_AF"),  
                                           valueName = "val", 
                                           argDesc = "The allele frequency cutoff to assign BA1 (benign) status.",
                                           defaultValue = Some(0.05)
                                           ) ::
                    new BinaryArgument[Double](name = "PM2_AF",
                                           arg = List("--PM2_AF"),  
                                           valueName = "val", 
                                           argDesc = "The allele frequency cutoff to assign PM2 (moderate pathogenic) status.",
                                           defaultValue = Some(0.0001)
                                           ) ::
                    //new BinaryOptionArgument[List[String]](
                    //                     name = "dropKeys", 
                    //                     arg = List("--dropKeys"), 
                    //                     valueName = "key1,key2,...",  
                    //                     argDesc =  "A list of INFO keys to omit from the output VCF."
                    //                    ) :: 
                    new FinalArgument[String](
                                         name = "invcf",
                                         valueName = "variants.vcf",
                                         argDesc = "input VCF file. This file must have been processed by " // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outvcf",
                                         valueName = "outvcf",
                                         argDesc = "The output vcf file."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
          
    def run(args : Array[String]) {
     val out = parser.parseArguments(args.toList.tail);
     if(out){
       
       CalcACMGVar.AssessACMGWalker(
                   chromList = parser.get[Option[List[String]]]("chromList"),
                   clinVarVcf = parser.get[Option[String]]("clinVarVcf").get,
                   txToGeneFile = parser.get[Option[String]]("txToGeneFile"),
                   ctrlAlleFreqKeys = parser.get[List[String]]("ctrlAlleFreqKeys"),
                   toleranceFile = parser.get[Option[String]]("toleranceFile"),
                   domainFile = parser.get[Option[String]]("domainFile"),
                   conservedElementFile = parser.get[Option[String]]("conservedElementFile"),
                   BA1_AF = parser.get[Double]("BA1_AF"),
                   PM2_AF = parser.get[Double]("PM2_AF"),
                   rmskFile = parser.get[Option[String]]("rmskFile")//,
                   //dropKeys = parser.get[Option[List[String]]]("dropKeys"),
                   ).walkVCFFile(
                   infile    = parser.get[String]("invcf"),
                   outfile   = parser.get[String]("outvcf"),
                   chromList = parser.get[Option[List[String]]]("chromList")
                   )
           
     }   
  }
  }
  
  case class AssessACMGWalker(
                 chromList : Option[List[String]],
                 clinVarVcf : String,
                 txToGeneFile : Option[String],
                 ctrlAlleFreqKeys : List[String],
                 toleranceFile : Option[String],
                 domainFile : Option[String],
                 conservedElementFile : Option[String],
                 BA1_AF : Double,
                 PM2_AF : Double,
                 rmskFile : Option[String]
                ) extends internalUtils.VcfTool.VCFWalker {
    
    reportln("Creating AssessACMGWalker() ["+stdUtils.getDateAndTimeString+"]","note")
    
    val inSilicoKeys : Seq[(String,Set[String],Set[String])] = Vector[(String,Set[String],Set[String])](
                        ("dbNSFP_MetaSVM_pred",Set[String]("D"),Set[String]("T"))//, //Includes predictors:   PolyPhen-2, SIFT, LRT, MutationTaster, Mutation Assessor, FATHMM, GERP++, PhyloP, SiPhy
                        //("dbNSFP_PROVEAN_pred",Set[String]("D"))
                    )
    val inSilicoMinCt = -1;
    val inSilicoMin = if(inSilicoMinCt == -1) inSilicoKeys.size else inSilicoMinCt;
    
    val txToGene : (String => String) = txToGeneFile match {
      case Some(f) => {
        val txToGeneMap = getLinesSmartUnzip(f).map(line => {
          val cells = line.split("\t");
          (cells(0),cells(1));
        }).toMap
        ((s : String) => txToGeneMap.getOrElse(s,s));
      }
      case None => {
        ((s : String) => s);
      }
    }

    
//getClinVarVariants(infile : String, chromList : Option[List[String]], vcfCodes : VCFAnnoCodes = VCFAnnoCodes()) : scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo)]]
    /****************************************
     * GENE MAPS:
     */ 
    reportln("Reading gene maps... ["+stdUtils.getDateAndTimeString+"]","debug");
    val geneIsLofSensitive : (String => Boolean) = toleranceFile match {
      case Some(f) => {
        val lines = getLinesSmartUnzip(f);
        val headerCells = lines.next.split("\\s+");
        val geneCol = headerCells.indexWhere(headerCell => headerCell == "geneID");
        val lofCol = headerCells.indexWhere(headerCell => headerCell == "LOFtolerant");
        var lofTolerantGeneList = Set[String]();
        while(lines.hasNext){
          val cells = lines.next().split("\\s+");
          if(cells(lofCol) == "1"){
            lofTolerantGeneList = lofTolerantGeneList + cells(geneCol);
          }
        }
        ((s : String) => {
          ! lofTolerantGeneList.contains(s);
        })
      }
      case None => {
        ((s : String) => true)        
      }
    }
    
    val geneIsMisSensitive : (String => Boolean) = toleranceFile match {
      case Some(f) => {
        val lines = getLinesSmartUnzip(f);
        val headerCells = lines.next.split("\\s+");
        val geneCol = headerCells.indexWhere(headerCell => headerCell == "geneID");
        val tolCol = headerCells.indexWhere(headerCell => headerCell == "LOFtolerant");
        var tolerantGeneList = Set[String]();
        while(lines.hasNext){
          val cells = lines.next().split("\\s+");
          if(cells(tolCol) == "1"){
            tolerantGeneList = tolerantGeneList + cells(geneCol);
          }
        }
        ((s : String) => {
          ! tolerantGeneList.contains(s);
        })
      }
      case None => {
        ((s : String) => true)        
      }
    }
    
    val geneIsMisInsensitive : (String => Boolean) = ((s : String) => { geneIsLofSensitive(s); })
    reportln("done with gene maps...["+stdUtils.getDateAndTimeString+"]","debug");
    //val geneIsMisSensitive : (String => Boolean) = ((s : String) => true)
    //val geneIsMisInsensitive : (String => Boolean) = ((s : String) => false)
    
    /****************************************
     * LOCUS MAPS:
     */
    val conservedArray = conservedElementFile match {
      case Some(f) => {
        reportln("reading conserved locus file ["+stdUtils.getDateAndTimeString+"]","debug");
        val arr : genomicAnnoUtils.GenomicArrayOfSets[String] = genomicAnnoUtils.GenomicArrayOfSets[String](false);
        val lines = getLinesSmartUnzip(f);
        
        lines.foreach(line => {
          val cells = line.split("\t");
          val (chrom,start,end) = (cells(0),string2int(cells(1)),string2int(cells(2)))
          arr.addSpan(commonSeqUtils.GenomicInterval(chrom, '.', start,end), "CE");
        })
        arr.finalizeStepVectors;
        reportln("done with conserved locus file ["+stdUtils.getDateAndTimeString+"]","debug");
        Some(arr);
      }
      case None => {
        None;
      }
    }
    val locusIsConserved  : (commonSeqUtils.GenomicInterval => Boolean) = conservedArray match {
      case Some(arr) => {
        (iv : commonSeqUtils.GenomicInterval) => {
          ! arr.findIntersectingSteps(iv).foldLeft(Set[String]()){case (soFar,(iv,currSet)) => {
            soFar ++ currSet;
          }}.isEmpty
        }
      }
      case None => {
        (iv : commonSeqUtils.GenomicInterval) => {
          true;
        }
      }
    }
    
    val locusArray = domainFile match {
      case Some(f) => {
        reportln("reading domain locus file ["+stdUtils.getDateAndTimeString+"]","debug");
        val lines = getLinesSmartUnzip(f);
        val headerCells = lines.next.split("\t");
        val chromCol = 0;
        val startCol = headerCells.indexWhere(headerCell => headerCell == "start");
        val endCol = headerCells.indexWhere(headerCell => headerCell == "end");
        val domainIdCol = headerCells.indexWhere(headerCell => headerCell == "domainID");
        
        val geneArray : genomicAnnoUtils.GenomicArrayOfSets[String] = genomicAnnoUtils.GenomicArrayOfSets[String](false);
        
        while(lines.hasNext){
          val cells = lines.next.split("\t");
          val (chrom,start,end,domainID) = (cells(chromCol),string2int(cells(startCol)),string2int(cells(endCol)),cells(domainIdCol))
          geneArray.addSpan(commonSeqUtils.GenomicInterval(chrom, '.', start,end), domainID);
        }
        geneArray.finalizeStepVectors;
        reportln("done with domain locus file. ["+stdUtils.getDateAndTimeString+"]","debug");
        Some(geneArray);
      }
      case None => {
        None;
      }
    }
    
    val locusDomains : ((commonSeqUtils.GenomicInterval) => Set[String]) = locusArray match {
      case Some(arr) => {
        (iv : commonSeqUtils.GenomicInterval) => {
          arr.findIntersectingSteps(iv).foldLeft(Set[String]()){case (soFar,(iv,currSet)) => {
            soFar ++ currSet;
          }}
        }
      }
      case None => {
        (iv : commonSeqUtils.GenomicInterval) => {
          Set[String]("UNK")
        }
      }
    }
    
    val locusIsHotspot    : (commonSeqUtils.GenomicInterval => Boolean) = (iv : commonSeqUtils.GenomicInterval) => {
      ! locusDomains(iv).isEmpty;
    }  //((s : String,i : Int) => false)

    
    val repLocusArray = rmskFile match {
      case Some(f) => {
        reportln("reading rmsk file ["+stdUtils.getDateAndTimeString+"]","debug");
        val geneArray : genomicAnnoUtils.GenomicArrayOfSets[String] = genomicAnnoUtils.GenomicArrayOfSets[String](false);
        getLinesSmartUnzip(f).foreach(line => {
          val cells = line.split("\t");
          val chrom = cells(5);
          val start = string2int(cells(6));
          val end = string2int(cells(7));
          geneArray.addSpan(commonSeqUtils.GenomicInterval(chrom,'.',start,end),cells(11));
        })
        geneArray.finalizeStepVectors;
        reportln("done with rmsk file ["+stdUtils.getDateAndTimeString+"]","debug");
        Some(geneArray);
      }
      case None => {
        None;
      }
    }
    
    val locusIsRepetitive : (commonSeqUtils.GenomicInterval => Boolean) = repLocusArray match {
      case Some(arr) => {
        ((iv : commonSeqUtils.GenomicInterval) => {
          arr.findIntersectingSteps(iv).foldLeft(Set[String]()){case (soFar,(iv,currSet)) => { soFar ++ currSet }}.isEmpty;
        })
      }
      case None => {
        ((iv : commonSeqUtils.GenomicInterval) => false)
      }
    }
    
    /****************************************
     * WALK FCN:
     */
    
    var vcfCodes : VCFAnnoCodes = VCFAnnoCodes();
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      val (clinVarVariantSet,clinVarVariants) : (scala.collection.Map[String,String],scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]]) = getClinVarVariants(infile=clinVarVcf,chromList =chromList, vcfCodes = vcfCodes);

      val newHeaderLines = List(
            new VCFInfoHeaderLine(vcfCodes.assess_IsRepetitive, 1, VCFHeaderLineType.Integer,    "Indicates whether the variant intersects with a repetitive region, as defined by repeatmasker (derived from rmsk.txt, downloaded from ucsc, APR 2017)"),
            new VCFInfoHeaderLine(vcfCodes.assess_IsConserved, 1, VCFHeaderLineType.Integer,    "Indicates whether the variant intersects with a disproportionately-conserved region (as defined by GERPplusplus)."),
            new VCFInfoHeaderLine(vcfCodes.assess_IsHotspot, 1, VCFHeaderLineType.Integer,    "Indicates whether the variant intersects with a known domain."),
            //new VCFInfoHeaderLine(vcfCodes.assess_domain, 1, VCFHeaderLineType.Integer,    ""),
            //new VCFInfoHeaderLine(vcfCodes.assess_criteria, 1, VCFHeaderLineType.Integer,    ""),
            //new VCFInfoHeaderLine(vcfCodes.assess_GeneSenseLOF, 1, VCFHeaderLineType.Integer,    ""),
            //new VCFInfoHeaderLine(vcfCodes.assess_GeneSenseMis, 1, VCFHeaderLineType.Integer,    ""),
            //new VCFInfoHeaderLine(vcfCodes.assess_ctrlAFMIN, 1, VCFHeaderLineType.Float,    ""),
            new VCFInfoHeaderLine(vcfCodes.assess_ctrlAFMAX, 1, VCFHeaderLineType.Float,    "The highest alt allele frequency found across the control datasets ("+ctrlAlleFreqKeys.mkString(",")+")"),
            new VCFInfoHeaderLine(vcfCodes.assess_LofGenes, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "A list of genes for which this variant causes loss-of-function (fs, early-stop, start-loss)"),
            new VCFInfoHeaderLine(vcfCodes.assess_MisGenes, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "A list of genes for which this variant causes a missense."),
            
            new VCFInfoHeaderLine(vcfCodes.assess_PVS1, 1, VCFHeaderLineType.Integer,    "Loss-of-function variant in gene that is LOF-sensitive. Loss-of-function is defined as one of:"+
                                                                                         "stop-gain, frameshift, total-gene-indel, or splice junction indel. "+
                                                                                         "A gene is defined as LOF-sensitive if at least 10pct of pathogenic "+
                                                                                         "clinvar variants are LOF-type variants."),
            new VCFInfoHeaderLine(vcfCodes.assess_PS1, 1, VCFHeaderLineType.Integer,      "Variant has the same amino acid change as a pathogenic variant from ClinVar."),
            new VCFInfoHeaderLine(vcfCodes.assess_PM1, 1, VCFHeaderLineType.Integer,      "Located in a known domain. (currently just any domain)"),
            new VCFInfoHeaderLine(vcfCodes.assess_PM2, 1, VCFHeaderLineType.Integer,      "Alt allele frequency less than or equal to "+PM2_AF+" in all control datasets ("+ctrlAlleFreqKeys.mkString(",")+")"),
            new VCFInfoHeaderLine(vcfCodes.assess_PM4, 1, VCFHeaderLineType.Integer,      "Protein length change in nonrepeat region OR stop-loss variant."),
            new VCFInfoHeaderLine(vcfCodes.assess_PM5, 1, VCFHeaderLineType.Integer,      "Novel missense variant change at amino acid where a different amino acid change is known pathogenic in ClinVar."),
            new VCFInfoHeaderLine(vcfCodes.assess_PP2, 1, VCFHeaderLineType.Integer,      "Missense variant in gene that is missense-sensitive (missense variants are at least 10pct of known pathogenic variants in clinvar)."),
            new VCFInfoHeaderLine(vcfCodes.assess_PP3, 1, VCFHeaderLineType.Integer,      "Predicted to be damaging by at least "+inSilicoMin+" of the following in silico prediction algorithms: "+inSilicoKeys.map(_._1).mkString(",")),
            new VCFInfoHeaderLine(vcfCodes.assess_BP1, 1, VCFHeaderLineType.Integer,      "Missense variant in gene that is NOT missense-sensitive (less than 10pct of pathogenic variants are missense)"),
            new VCFInfoHeaderLine(vcfCodes.assess_BP3, 1, VCFHeaderLineType.Integer,      "In-frame indels in a repetitive region that does NOT overlap with any known domain."),
            new VCFInfoHeaderLine(vcfCodes.assess_BP4, 1, VCFHeaderLineType.Integer,      "Predicted to be benign by at least "+inSilicoMin+" of the following in silico prediction algorithms: "+inSilicoKeys.map(_._1).mkString(",")),
            new VCFInfoHeaderLine(vcfCodes.assess_BP7, 1, VCFHeaderLineType.Integer,      "Synonymous variant that does NOT intersect with a conserved element region."),
            new VCFInfoHeaderLine(vcfCodes.assess_BA1,    1, VCFHeaderLineType.Integer,   "Allele frequency greater than 5 percent in one or more control dataset ("+ctrlAlleFreqKeys.mkString(",")+")"),
            new VCFInfoHeaderLine(vcfCodes.assess_RATING, 1, VCFHeaderLineType.String,    "ACMG Pathogenicity rating: PATHO - pathogenic. LPATH - likely pathogenic, VUS - variant, unknown significance, LB - likely benign, B - benign."),
            new VCFInfoHeaderLine(vcfCodes.assess_WARNFLAG, 1, VCFHeaderLineType.Integer, "Whether or not there is anything odd about this variant that may require manual inspection."),
            new VCFInfoHeaderLine(vcfCodes.assess_WARNINGS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "List of warnings concerning this variant."),
            
            new VCFInfoHeaderLine(vcfCodes.geneIDs, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "List of gene symbols that this variant overlaps with."),
            
            new VCFInfoHeaderLine(vcfCodes.assess_exactMatchRS, 1, VCFHeaderLineType.String,   "RS number referring to pathogenic ClinVar variant that matches this variant exactly (or blank if there is no matching variant)."),
            new VCFInfoHeaderLine(vcfCodes.assess_aminoMatchRS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "List of RS numbers referring to pathogenic ClinVar variant(s) that match the amino acid change of this variant  (or blank if there is no matching variant)."),
            new VCFInfoHeaderLine(vcfCodes.assess_nearMatchRS,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "List of RS numbers referring to pathogenic ClinVar variant(s) that alter the same amino acid position, but change to a different amino acid  (or blank if there is no matching variant).")
          );
            
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      reportln("Walking input VCF...","note")
      return ( 
      vcIter.map(v => {
        assessVariant(v=v,
            clinVarVariants=clinVarVariants,
            txToGene = txToGene,
            geneIsLofSensitive = geneIsLofSensitive,
            geneIsMisSensitive = geneIsMisSensitive,
            geneIsMisInsensitive = geneIsMisInsensitive,
            locusIsRepetitive = locusIsRepetitive,
            locusIsConserved = locusIsConserved,
            locusIsHotspot = locusIsHotspot,
            ctrlAlleFreqKeys = ctrlAlleFreqKeys,
            inSilicoKeys = inSilicoKeys,
            inSilicoMinCt = inSilicoMinCt,
            BA1_AF = BA1_AF,
            PM2_AF = PM2_AF,
            clinVarVariantSet = clinVarVariantSet,
            vcfCodes = vcfCodes)
       }), newHeader );
    }
    reportln("AssessACMGWalker() Created... ["+stdUtils.getDateAndTimeString+"]","note")
  }
  
  def assessVariant(v : VariantContext, 
                    clinVarVariants : scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]],
                    txToGene : (String => String) = ((s : String) => s),
                    geneIsLofSensitive : (String => Boolean) = ((s : String) => true),
                    geneIsMisSensitive : (String => Boolean) = ((s : String) => true),
                    geneIsMisInsensitive : (String => Boolean) = ((s : String) => false),
                    locusIsRepetitive : (commonSeqUtils.GenomicInterval => Boolean) ,
                    locusIsConserved  : (commonSeqUtils.GenomicInterval => Boolean) ,
                    locusIsHotspot    : (commonSeqUtils.GenomicInterval => Boolean) ,
                    ctrlAlleFreqKeys : Seq[String] = Seq("1KG_AF","ESP_EA_AF","ExAC_ALL","SWH_AF_GRP_CTRL"),
                    
                    inSilicoKeys : Seq[(String,Set[String],Set[String])] = Vector[(String,Set[String],Set[String])](
                        ("dbNSFP_MetaSVM_pred",Set[String]("D"),Set[String]("T"))//, //Includes predictors:   PolyPhen-2, SIFT, LRT, MutationTaster, Mutation Assessor, FATHMM, GERP++, PhyloP, SiPhy
                        //("dbNSFP_PROVEAN_pred",Set[String]("D"))
                    ),
                    inSilicoMinCt : Int = -1, //-1 == ALL
                    inSilicoToleratedCt : Int = 1,
                    
                    BA1_AF : Double = 0.05,
                    PM2_AF : Double = 0.0001,
                    clinVarVariantSet : scala.collection.Map[String,String],
                    vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                    
                    //clinVarVariants?
                    ) : VariantContext = {
    
    var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
    val crits = new ACMGCritSet();

    var problemList = Set[String]();
    
    val inSilicoMin = if(inSilicoMinCt == -1) inSilicoKeys.size else inSilicoMinCt;
    
    val chrom = v.getContig();
    val pos = v.getStart();
    
    val refAlle = v.getReference();
    val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a)).zipWithIndex.filter{case (alt,altIdx) => { alt.getBaseString() != "*" }}

    val txList = v.getAttributeAsList(vcfCodes.txList_TAG).toVector.map(_.toString).filter(_ != ".");
    
    if(txList.length > 0){
    val geneList = txList.map(x => txToGene(x));
    vb = vb.attribute(vcfCodes.geneIDs , geneList.mkString(","));
    
    val vTypesList = v.getAttributeAsList(vcfCodes.vType_TAG).toVector.map(_.toString.split("\\|").toVector);
    val vLVL       = v.getAttributeAsList(vcfCodes.vMutLVL_TAG).toVector.map(_.toString);
    val vMutPList = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.map(_.toString.split("\\|").toVector);
    val vMutCList = v.getAttributeAsList(vcfCodes.vMutC_TAG).toVector.map(_.toString.split("\\|").toVector);
    val vMutInfoList = v.getAttributeAsList(vcfCodes.vMutINFO_TAG).toVector.map(_.toString.split("\\|").toVector.map(x => {internalUtils.TXUtil.getPvarInfoFromString(x)}));
    
    val variantIV = commonSeqUtils.GenomicInterval(v.getContig(),'.', start = v.getStart() - 1, end = math.max(v.getEnd(),v.getStart+1));
    
    val isRepetitive = locusIsRepetitive(variantIV);
    val isConserved  = locusIsConserved(variantIV);
    val isHotspot    = locusIsHotspot(variantIV);
    
    if(altAlleles.length > 1){
      error("FATAL ERROR: Cannot deal with multiallelic variants! Split variants to single-allelic!");
    }
    
    vb = vb.attribute(vcfCodes.assess_IsRepetitive , if(isRepetitive) "1" else "0");
    vb = vb.attribute(vcfCodes.assess_IsConserved ,  if(isConserved)  "1" else "0");
    vb = vb.attribute(vcfCodes.assess_IsHotspot ,    if(isHotspot)    "1" else "0");
    //vb = vb.attribute(vcfCodes.assess_domain ,       if(isHotspot)    "1" else "0");
    
    val (alt,altIdx) = altAlleles.head;
    
    //val out : Vector[ACMGCritSet] = altAlleles.map{ case (alt,altIdx) => {
      val vTypes = vTypesList(altIdx).map(_.split("_"));
      val vMutP = vMutPList(altIdx);
      val vMutC = vMutCList(altIdx);
      val vInfo = vMutInfoList(altIdx);
      val combo = txList.zip(vInfo).zip(vMutC).zipWithIndex.map{case (((tx,info),c),i) => (txToGene(tx),tx,info,c,i)}
      
      if(combo.exists{ case (g,tx,info,c,i) => { info.subType == "???" }}){
        problemList = problemList + "VT|OddballVariantType" + "VT|???";
      }
      if(combo.exists{ case (g,tx,info,c,i) => { info.subType == "FullExonIndel" }}){
        problemList = problemList + "VT|OddballVariantType" + "VT|FullExonIndel";
      }
      if(combo.exists{ case (g,tx,info,c,i) => { info.subType == "FullIntronIndel" }}){
        problemList = problemList + "VT|OddballVariantType" + "VT|FullIntronIndel";
      }
      if(combo.exists{ case (g,tx,info,c,i) => { info.subType == "total-loss" }}){
        problemList = problemList + "VT|OddballVariantType" + "VT|totalloss";
      }
      
      //***************************HotSpot:
      vb =  vb.attribute(vcfCodes.assess_PM1 , if(isHotspot) "1" else "0");
      crits.addCrit(new ACMGCrit("PM",1,isHotspot,Seq[String]()))
      
      //*****************ALLELE FREQS:
      val splitIdx = v.getAttributeAsInt(vcfCodes.splitIdx_TAG,0);
      val numSplit = v.getAttributeAsInt(vcfCodes.numSplit_TAG,1);
      val ctrlAlleFreqs = ctrlAlleFreqKeys.map(key => {
        val afList = v.getAttributeAsList(key).map(x => x.toString()).map(x => {
          if(x == ".") 0.toDouble;
          else string2double(x);
        })
        if(afList.length == 0){
          0.toDouble;
        } else if(afList.length == 1){
          afList(0);
        } else if(afList.length == numSplit){
          afList(splitIdx);
        } else {
          problemList = problemList + "AF|popAFalleleMismatch";
          warning("Warning: allele freq annotation doesn't match numSplit: afList.length = "+afList.length+", numSplit = "+numSplit+"\n"+
                  "   "+"ref["+refAlle.getBaseString()+"],ALTs=["+altAlleles.map(_._1.getBaseString()).mkString(",")+"]"+key+"=["+v.getAttributeAsList(key).map(_.toString()).mkString(",")+"]\n"+
                  "   "+v.toStringWithoutGenotypes(),"POPAF_FORMAT_WARNING",25);
          0.toDouble;
        }
      })

      //val (minAF, maxAF) = ctrlAlleFreqs.foldLeft( (1.toDouble,0.toDouble) ){ case ((minSF,maxSF),af) => {
      //  
       // val afOpt = if(afString == "."){
      //    None
       // } else Some(string2double(afString));
      //  (minOption(Some(minSF),afOpt).get, maxOption(Some(maxSF),afOpt).get)
      //}}
      val maxAF = ctrlAlleFreqs.foldLeft(0.toDouble){case (maxSF,af) => {
        math.max(maxSF,af);
      }}
      //vb =  vb.attribute(vcfCodes.assess_ctrlAFMIN , minAF.toString);
      vb =  vb.attribute(vcfCodes.assess_ctrlAFMAX , maxAF.toString);
      vb =  vb.attribute(vcfCodes.assess_BA1 , if(maxAF >= BA1_AF) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_PM2 , if(maxAF <= PM2_AF) "1" else "0");
      crits.addCrit(new ACMGCrit("BA",1,maxAF >= BA1_AF,Seq[String]()));
      crits.addCrit(new ACMGCrit("PM",2,maxAF <= PM2_AF,Seq[String]()));
      
      //***************************LOF:
      val LofGenes = combo.filter{case (g,tx,info,c,i) => {
        info.severityType == "LLOF" || info.subType == "START-LOSS" || info.subType == "splice"
      }}.map{case (g,tx,info,c,i) => {
        g
      }}.toSet.toList.sorted;
      vb =  vb.attribute(vcfCodes.assess_LofGenes , LofGenes.mkString(","));
      val pvs1flag = LofGenes.exists(g => {
        geneIsLofSensitive(g);
      })
      vb =  vb.attribute(vcfCodes.assess_PVS1 , if(pvs1flag) "1" else "0");
      crits.addCrit(new ACMGCrit("PVS",1,pvs1flag,Seq[String]()));
      //***************************Mis:
      val MisGenes = combo.filter{case (g,tx,info,c,i) => {
        info.subType == "swapAA";
      }}.map{case (g,tx,info,c,i) => {
        g
      }}.toSet.toList.sorted;
      vb =  vb.attribute(vcfCodes.assess_MisGenes , MisGenes.mkString(","));
      val pp2flag = MisGenes.exists(g => {
        geneIsMisSensitive(g);
      });
      val bp1flag = (! pvs1flag) && (MisGenes.length > 0) && MisGenes.forall(g => {
        geneIsMisInsensitive(g);
      })
      vb =  vb.attribute(vcfCodes.assess_PP2 , if(pp2flag) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_BP1 , if(bp1flag) "1" else "0");
      crits.addCrit(new ACMGCrit("PP",2,pp2flag,Seq[String]()));
      crits.addCrit(new ACMGCrit("BP",1,bp1flag,Seq[String]()));
      //******************************* Length change variants:
      val pm4flag = {
        ((! isRepetitive) && ({
           combo.exists{case (g,tx,info,c,i) => {
             info.subType == "insAA" || info.subType == "delAA" || info.subType == "indelAA"
           }}
        })) || combo.exists{case (g,tx,info,c,i) => {
             info.pType == "STOP-LOSS"
           }}
      }
      vb =  vb.attribute(vcfCodes.assess_PM4,if(pm4flag) "1" else "0");
      val bp3flag = {
        (isRepetitive) && (! isHotspot) && ({
           combo.exists{case (g,tx,info,c,i) => {
             info.subType == "insAA" || info.subType == "delAA" || info.subType == "indelAA"
           }}
        })
      }
      vb =  vb.attribute(vcfCodes.assess_BP3,if(bp3flag) "1" else "0");
      
      crits.addCrit(new ACMGCrit("BP",3,bp3flag,Seq[String]()));
      crits.addCrit(new ACMGCrit("PM",4,pm4flag,Seq[String]()));
      
      //******************************* ClinVar matching:
      
      var ps1 = ACMGCrit("PS",1,false,Seq[String]());
      var pm5 = ACMGCrit("PM",5,false,Seq[String]());
      var ps1RS = Set[String]();
      var pm5RS = Set[String]();
      val variantString = v.getContig() + ":"+v.getAttributeAsString(vcfCodes.vMutG_TAG,"?!?");
      val exactRS = clinVarVariantSet.getOrElse(variantString,"")
      
      combo.filter{case (g,tx,info,c,i) => {
        info.severityType == "NONSYNON"
      }}.foreach{case (g,tx,info,c,i) => {
        val cvExactMatch = clinVarVariants(tx).exists{case (cvc,cvinfo,rsnum) => {cvinfo.start == info.start && cvinfo.subType == info.subType && cvinfo.pvar == info.pvar}};
        if(cvExactMatch){
          ps1 = ps1.merge(ACMGCrit("PS",1,true,Seq[String](tx + ":" + info.pvar)));
          ps1RS = ps1RS ++ clinVarVariants(tx).filter{ case (cvc,cvinfo,rsnum) => {cvinfo.start == info.start && cvinfo.subType == info.subType && cvinfo.pvar == info.pvar}}.map(_._3).toSet;
        } else if(info.subType == "swapAA"){
            val partialMatch = clinVarVariants(tx).filter{case (cvc,cvinfo,rsnum) => {cvinfo.subType == "swapAA" && cvinfo.start == info.start && cvinfo.altAA != info.altAA}};
            if(partialMatch.size > 0){
              pm5 = pm5.merge(ACMGCrit("PM",5,true,Seq[String](tx+":"+info.pvar+"Vs"+partialMatch.head._2.altAA)));
              pm5RS = pm5RS ++ partialMatch.map(_._3).toSet;
            }
        }
      }}
      if(ps1.flag){
        vb =  vb.attribute(vcfCodes.assess_PS1,"1");
        vb =  vb.attribute(vcfCodes.assess_PM5,"0");
      } else if(pm5.flag){
        vb =  vb.attribute(vcfCodes.assess_PS1,"0");
        vb =  vb.attribute(vcfCodes.assess_PM5,"1");
      } else {
        vb =  vb.attribute(vcfCodes.assess_PS1,"0");
        vb =  vb.attribute(vcfCodes.assess_PM5,"0");
      }
      crits.addCrit(ps1);
      crits.addCrit(pm5);
      
      vb = vb.attribute(vcfCodes.assess_exactMatchRS, exactRS);
      vb = vb.attribute(vcfCodes.assess_aminoMatchRS, ps1RS.toVector.sorted.mkString(","));
      vb = vb.attribute(vcfCodes.assess_nearMatchRS,  pm5RS.toVector.sorted.mkString(","));
      
      //******************************* Insilico:
      //inSilicoMin
      //val numSplit = v.getAttributeAsInt(vcfCodes.numSplit_TAG,1);
      val rawAlles = v.getAttributeAsList(vcfCodes.splitAlle_TAG).map(_.toString());
      
      val isBadFlag = if(alt.getBaseString().length > 1){
        false;
      } else if(numSplit == 1){
        var warn : Seq[String] = Seq[String]();
        inSilicoKeys.count{case (key,damSet,safeSet) => {
          val xraw = v.getAttributeAsList(key);
          if(xraw.length == 1){
            damSet.contains(xraw.head.toString());
          } else if(xraw.length == 0){
            false;
          } else {
            warning("One known allele, multiple InSilico fields for key: "+key+"\n"+"line:\n   "+v.toString(), "maybeAmbigInSilicoWarning", limit = 25);
            problemList = problemList + "IS|inSilicoAlleleMismatch1"
            xraw.map(_.toString()).exists(damSet.contains(_));
          }
        }} >= inSilicoMin;
      } else {
        val inSilicoRaw = inSilicoKeys.map{case (key,damSet,safeSet) => {
          v.getAttributeAsList(key).toList;
        }}
        if(rawAlles.length != numSplit){
          //warning("", "", limit = 25);
          error("Incompatible VCF INFO values for tags "+vcfCodes.numSplit_TAG +" and " +vcfCodes.splitAlle_TAG+"!\n"+
                "numSplit = "+numSplit+"\n"+
                "rawAlles = "+rawAlles.mkString(",")+"\n"
          );
        }
        if(splitIdx >= numSplit){
          error("Incompatible VCF INFO values for tags "+vcfCodes.numSplit_TAG +" and " +vcfCodes.splitIdx_TAG+"!\n"+
                "numSplit = "+numSplit+"\n"+
                "splitIdx = "+splitIdx+"\n"
          );
        }
        val keepAlles = rawAlles.zipWithIndex.filter{case (a,i) =>{
          a.length == 1
        }}
        keepAlles.zipWithIndex.find{case ((a,oldidx),newidx) => {
          oldidx == splitIdx;
        }} match {
          case Some(((a,oldidx),newidx)) => {
            inSilicoKeys.count{case (key,damSet,safeSet) => {
                val xraw = v.getAttributeAsList(key);
                if(xraw.length == keepAlles.length){
                  damSet.contains(xraw(newidx).toString());
                } else if(xraw.length == 0){
                  false;
                } else {
                  problemList = problemList + "IS|inSilicoAlleleMismatch2"
                  warning("Multiple known alleles, wrong # of InSilico fields for key: "+key+"\n"+""+
                          "Alt alleles (in order) = "+Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a)).map(_.getBaseString()).mkString(",")+", rawAlles="+rawAlles.mkString(",")+", InSilico("+key+")="+xraw.map(_.toString()).mkString(",")+"\n"+
                          "line:\n"+v.toString(), "ambigInSilicoWarning", limit = 100);
                  false;
                }
            }} >= inSilicoMin;
          }
          case None => {
            false;
          }
        }
      }
      val isOkFlag = if(alt.getBaseString().length > 1){
        false;
      } else if(numSplit == 1){
        var warn : Seq[String] = Seq[String]();
        inSilicoKeys.count{case (key,damSet,safeSet) => {
          val xraw = v.getAttributeAsList(key);
          if(xraw.length == 1){
            safeSet.contains(xraw.head.toString());
          } else if(xraw.length == 0){
            false;
          } else {
             problemList = problemList + "IS|inSilicoAlleleMismatch1"
             warning("One known allele, wrong # of InSilico fields for key: "+key+"\n"+""+
                     "Alt alleles (in order) = "+Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a)).map(_.getBaseString()).mkString(",")+", rawAlles="+rawAlles.mkString(",")+", InSilico("+key+")="+xraw.map(_.toString()).mkString(",")+"\n"+
                     "line:\n"+v.toString(), "maybeAmbigInSilicoWarning", limit = 100);
            //warning("One known allele, multiple InSilico fields for key: "+key+"\n"+"line:\n   "+v.toString(), "maybeAmbigInSilicoWarning", limit = 25);
            xraw.map(_.toString()).exists(safeSet.contains(_));
          }
        }} >= inSilicoToleratedCt;
      } else {
        val inSilicoRaw = inSilicoKeys.map{case (key,damSet,safeSet) => {
          v.getAttributeAsList(key).toList;
        }}
        if(rawAlles.length != numSplit){
          //warning("", "", limit = 25);
          error("Incompatible VCF INFO values for tags "+vcfCodes.numSplit_TAG +" and " +vcfCodes.splitAlle_TAG+"!\n"+
                "numSplit = "+numSplit+"\n"+
                "rawAlles = "+rawAlles.mkString(",")+"\n"
          );
        }
        if(splitIdx >= numSplit){
          error("Incompatible VCF INFO values for tags "+vcfCodes.numSplit_TAG +" and " +vcfCodes.splitIdx_TAG+"!\n"+
                "numSplit = "+numSplit+"\n"+
                "splitIdx = "+splitIdx+"\n"
          );
        }
        val keepAlles = rawAlles.zipWithIndex.filter{case (a,i) =>{
          a.length == 1
        }}
        keepAlles.zipWithIndex.find{case ((a,oldidx),newidx) => {
          oldidx == splitIdx;
        }} match {
          case Some(((a,oldidx),newidx)) => {
            inSilicoKeys.count{case (key,damSet,safeSet) => {
                val xraw = v.getAttributeAsList(key);
                if(xraw.length == keepAlles.length){
                  safeSet.contains(xraw(newidx).toString());
                } else if(xraw.length == 0){
                  false;
                } else {
                  problemList = problemList + "IS|inSilicoAlleleMismatch2"
                  warning("Multiple known alleles, wrong # of InSilico fields for key: "+key+"\n"+"line:\n   "+v.toString(), "ambigInSilicoWarning", limit = 25);
                  false;
                }
            }} >= inSilicoToleratedCt;
          }
          case None => {
            false;
          }
        }
      }
      
      vb = vb.attribute(vcfCodes.assess_PP3, if(isBadFlag) "1" else "0" );
      vb = vb.attribute(vcfCodes.assess_BP4, if(isOkFlag) "1" else "0" );
      
      crits.addCrit(new ACMGCrit("PP",3,isBadFlag,Seq[String]()));
      crits.addCrit(new ACMGCrit("BP",4,isOkFlag,Seq[String]()));
      
      //******************************* BP7: Synonymous w/ no splice impact, not highly conserved:
      val bp7flag = (! isConserved) && combo.forall{case (g,tx,info,c,i) => {
        info.severityType == "SYNON" || info.severityType == "PSYNON" || info.severityType == "UNK"
      }}
      vb = vb.attribute(vcfCodes.assess_BP7, if(bp7flag) "1" else "0" );
      crits.addCrit(new ACMGCrit("BP",7,bp7flag,Seq[String]()));
      
      //*******************************
      
      for((critLvl, critNum, critTAG) <- VcfTool.UNAUTOMATED_ACMG_PARAMS){
        if(v.hasAttribute(critTAG)){
          crits.addCrit(new ACMGCrit(critLvl,critNum,v.getAttributeAsString(critTAG,"") == "1",Seq[String]()));
        }
      }
      
      val rating = CalcACMGVar.getACMGPathogenicityRating(crits.getCritSet);
      vb = vb.attribute(vcfCodes.assess_RATING, rating );
      
      
    //vb =  vb.attribute(vcfCodes.assess_PP3,if(inSilicoFlag,"1","0"));
      //val inSilico
      
      
      //BP3:
      
      //crits;
    //}}.toVector;
    
    //return (out,vb.make());
    }
    
    vb = vb.attribute(vcfCodes.assess_WARNFLAG, if(problemList.isEmpty) "0" else "1" );
    vb = vb.attribute(vcfCodes.assess_WARNINGS, problemList.toVector.sorted.mkString(",") );
    
    return vb.make();
  }
  

  //output: tx => Set[(HGVDc,pVariantInfo])]
  def getClinVarVariants(infile : String, chromList : Option[List[String]], vcfCodes : VCFAnnoCodes = VCFAnnoCodes()) : (scala.collection.Map[String,String],scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]]) = {
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(infile, 
                                                                  chromList = chromList,
                                                                  vcfCodes = vcfCodes);
    
    val out = new scala.collection.mutable.AnyRefMap[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]](((k : String) => Set[(String,internalUtils.TXUtil.pVariantInfo,String)]()));
    var gset = new scala.collection.mutable.AnyRefMap[String,String]();
    
    reportln("Starting VCF read/write...","progress");
    for(v <- vcIter){
      try {
      val refAlle = v.getReference();
      val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a));
      //val vTypesList = v.getAttributeAsList(vcfCodes.vType_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
      val txList = v.getAttributeAsList(vcfCodes.txList_TAG).toVector.map(_.toString).filter(_ != ".");
      //val vMutPList = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
      val vMutCListRaw = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.filter(_ != ".");
      val vMutCList = vMutCListRaw.map(_.toString.split("\\|").toVector);
      val chrom = v.getContig();
      
      val vMutGList = v.getAttributeAsList(vcfCodes.vMutG_TAG).toVector.filter(_ != ".");
      vMutGList.foreach(g => {
        gset(chrom + ":" + g) = v.getID();
      })
      
      val vMutInfoList = if(txList.length > 0) {
          v.getAttributeAsList(vcfCodes.vMutINFO_TAG).toVector.map( (attrObj) => {
          val attrString = attrObj.toString();
          val attrSplit = attrString.split("\\|").toVector
          //reportln("vMutInfoListDebug: attrString = \""+attrString+"\"","debug");
          //reportln("vMutInfoListDebug: attrSplit["+attrSplit.length+"] = [\""+attrSplit.mkString("\", \"")+"\"]","debug");
          
          attrSplit.map(x => {
            internalUtils.TXUtil.getPvarInfoFromString(x)
          }) 
        })
      } else {
        Vector();
      }
      val vClnSig = v.getAttributeAsList("CLNSIG").toVector.map(_.toString());
      
      
      if(txList.length > 0) { 
        for((alle,altIdx) <- altAlleles.zipWithIndex){
           //val vTypes = vTypesList(altIdx);
           //val vMutP = vMutPList(altIdx);
           val vMutC = vMutCList(altIdx);
           val vInfo = vMutInfoList(altIdx);
           val clnSig = vClnSig(altIdx).split("\\|");
           val hasBenign = clnSig.exists(p => p == "2" || p == "3");
           val hasPatho  = clnSig.exists(p => p == "4" || p == "5");
           
           if(vMutC.length != txList.length){
             reportln("vMutC.length = "+vMutC.length+", txList = "+txList+"\nvMutC = [\""+vMutC.mkString("\",\"")+"\"]","debug");
           }
           if(hasPatho && (! hasBenign)){
             for(i <- Range(0,txList.length)){
               val tx = txList(i);
               out(tx) = out(tx) + ((vMutC(i),vInfo(i),v.getID()));
             }
           }
        }
      }

      } catch {
        case e : Exception => {
          reportln("Caught Exception on line:","note");
          reportln(v.toStringWithoutGenotypes(),"note");
          throw e;
        }
      }
    }
    return (gset,out);
  }
  
  
  
//##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Variant names from HGVS.    The order of these variants corresponds to the order of the info in the other clinical  INFO tags.">
//##INFO=<ID=CLNALLE,Number=.,Type=Integer,Description="Variant alleles from REF or ALT columns.  0 is REF, 1 is the first ALT allele, etc.  This is used to match alleles with other corresponding clinic
//##INFO=<ID=CLNSRC,Number=.,Type=String,Description="Variant Clinical Chanels">
//##INFO=<ID=CLNORIGIN,Number=.,Type=String,Description="Allele Origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - matern
//##INFO=<ID=CLNSRCID,Number=.,Type=String,Description="Variant Clinical Channel IDs">
//##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6
//##INFO=<ID=CLNDSDB,Number=.,Type=String,Description="Variant disease database name">
//##INFO=<ID=CLNDSDBID,Number=.,Type=String,Description="Variant disease database ID">
//##INFO=<ID=CLNDBN,Number=.,Type=String,Description="Variant disease name">
//##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="no_assertion - No assertion provided, no_criteria - No assertion criteria provided, single - Criteria provided single submitter, mult - Criteria
//##INFO=<ID=CLNACC,Number=.,Type=String,Description="Variant Accession and Versions">
  
  case class ClinVariantHolder(txid : String, varp : String, varc : String, varType : Array[String]){
    lazy val isSwap : Boolean = varType.last == "swapAA";
    lazy val posAA : Int = if(isSwap){
      string2int(varc.slice(5,varc.length-3));
    } else {
      -1; //NOT IMPLEMENTED YET!
    }
  }
  
  def getRepetitiveFunction(rmskFile : String) : ((String,Int) => Boolean) = {
    val geneArray : internalUtils.genomicAnnoUtils.GenomicArrayOfSets[String] = internalUtils.genomicAnnoUtils.GenomicArrayOfSets[String](false);
    
    for(line <- getLinesSmartUnzip(rmskFile)){
      val cells = line.split("\t");
      val (chrom,start,end,repClass) = (cells(5),string2int(cells(6)),string2int(cells(7)),cells(11));
      geneArray.addSpan(internalUtils.commonSeqUtils.GenomicInterval(chrom,'.',start,end),repClass);
    }
    
    val finalArray = geneArray.finalizeStepVectors;
    
    return ((chrom : String, pos : Int) => {
      ! finalArray.findSetAtPosition(chrom,pos,'.').isEmpty
    })
  }
  
  def assessVariantOLD(v : VariantContext, 
                    clinVarVariants : scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo)]],
                    txIsLofSensitive : (String => Boolean) = ((s : String) => true),
                    txIsMisSensitive : (String => Boolean) = ((s : String) => true),
                    locusIsRepetitive : ((String,Int) => Boolean) = ((s : String,i : Int) => false),
                    locusIsConserved  : ((String,Int) => Boolean) = ((s : String,i : Int) => true ),
                    locusIsHotspot    : ((String,Int) => Boolean) = ((s : String,i : Int) => false ),
                    ctrlAlleFreqKeys : Seq[String] = Seq("1KG_AF","ESP_EA_AF","ExAC_ALL","SWH_AF_GRP_CTRL"),
                    
                    vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                    //clinVarVariants?
                    ) : (Vector[ACMGCritSet],VariantContext) = {
    var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
    
    val chrom = v.getContig();
    val pos = v.getStart();
    
    val refAlle = v.getReference();
    val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a)).zipWithIndex.filter{case (alt,altIdx) => { alt.getBaseString() != "*" }}
    
    val vTypesList = v.getAttributeAsList(vcfCodes.vType_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
    val vLVL       = v.getAttributeAsList(vcfCodes.vMutLVL_TAG).toVector.map(_.toString);
    val txList = v.getAttributeAsList(vcfCodes.txList_TAG).toVector.map(_.toString);
    val vMutPList = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
    val vMutCList = v.getAttributeAsList(vcfCodes.vMutC_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
    val vMutInfoList = v.getAttributeAsList(vcfCodes.vMutINFO_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector.map(x => {internalUtils.TXUtil.getPvarInfoFromString(x)}));
    
    val isRepetitive = locusIsRepetitive(v.getContig(),v.getStart());
    val isConserved  = locusIsConserved(v.getContig(),v.getStart());
    val isHotspot    = locusIsHotspot(v.getContig(),v.getStart());
    
    if(altAlleles.length > 1){
      error("FATAL ERROR: Cannot deal with multiallelic variants! Split variants to single-allelic!");
    }
    
    vb = vb.attribute(vcfCodes.assess_IsRepetitive , if(isRepetitive) "1" else "0");
    vb = vb.attribute(vcfCodes.assess_IsConserved ,  if(isConserved)  "1" else "0");
    vb = vb.attribute(vcfCodes.assess_IsHotspot ,    if(isHotspot)    "1" else "0");
    //vb = vb.attribute(vcfCodes.assess_domain ,       if(isHotspot)    "1" else "0");

    val (alt,altIdx) = altAlleles.head;
    
    //val out : Vector[ACMGCritSet] = altAlleles.map{ case (alt,altIdx) => {
      val crits = new ACMGCritSet();
      val vTypes = vTypesList(altIdx).map(_.split("_"));
      val vMutP = vMutPList(altIdx);
      val vMutC = vMutCList(altIdx);
      val vInfo = vMutInfoList(altIdx);
      
      val (pm2,ba1) = checkDbAlleFreq(v,altIdx);
      crits.addCrit(pm2);
      crits.addCrit(ba1);
      crits.addCrit(checkLOF(vTypes,txList,vMutP,txIsLofSensitive));
      
      var bp7 = ACMGCrit("BP",7,false,Seq[String]());
      if( vTypes.forall(vt => {vt(1) != "LLOF" && vt(1) != "PLOF" && vt(1) != "NONSYNON"}) ){
        if(locusIsConserved(chrom,pos)) bp7 = ACMGCrit("BP",7,true, Seq[String]());
        //else                              bp7 = ACMGCrit("BP",7,false,Seq[String]());
      }
      
      var pm4 = ACMGCrit("PM",4,false,Seq[String]());
      val pm4List = vTypes.filter(_(1) == "NONSYNON").filter(vt => vt(2) == "STOP-LOSS" || vt.last == "insAA" || vt.last == "delAA" || vt.last == "indelAA")
      if( pm4List.length > 0 ){
        if(! isRepetitive) pm4 = ACMGCrit("PM",4,true, Seq[String]());
      }
      val combo = txList.zip(vInfo).zip(vMutC).zipWithIndex.map{case (((tx,info),c),i) => (tx,info,c,i)}
      
      var ps1 = ACMGCrit("PS",1,false,Seq[String]());
      var pm5 = ACMGCrit("PM",5,false,Seq[String]());
      var bpa3 = ACMGCrit("BPa",3,false,Seq[String]());
      
      combo.foreach{case (tx,info,c,i) => {
        val cvExactMatch = clinVarVariants(tx).filter{case (cvc,cvinfo) => {cvinfo.start == info.start && cvinfo.subType == info.subType && cvinfo.pvar == info.pvar}};
        if(cvExactMatch.size > 0){
          ps1 = ps1.merge(ACMGCrit("PS",1,true,Seq[String](tx + ":" + info.pvar)));
        }
        if(info.subType == "swapAA"){
          if(cvExactMatch.size == 0){
            val partialMatch = clinVarVariants(tx).filter{case (cvc,cvinfo) => {cvinfo.subType == "swapAA" && cvinfo.start == info.start && cvinfo.altAA != info.altAA}};
            if(partialMatch.size > 0){
              pm5 = pm5.merge(ACMGCrit("PM",5,true,Seq[String](tx+":"+info.pvar+"Vs"+partialMatch.head._2.altAA)));
            }
          }
        }
        if( info.subType == "delAA" || info.subType == "insAA" || info.subType == "multSwapAA" || info.subType == "indelAA" ){
          if(isRepetitive) bpa3.merge(ACMGCrit("BP",3,true,Seq[String](tx+":"+info.pvar)));
        }
      }}

      crits.addCrit(ps1);
      crits.addCrit(pm5);
      crits.addCrit(bpa3);
      //BP3:
      
      //crits;
    //}}.toVector;
    
    //return (out,vb.make());
      
    return null;
  }
  
  //def checkInFrameIndels
  
  def checkLOF(vTypes : Vector[Array[String]], txList : Vector[String], vMutP : Vector[String], txIsLofSensitive : (String => Boolean) = ((s : String) => true)) : ACMGCrit = {
    vTypes.zip(txList).zip(vMutP).map{case ((vtype,tx),aa) => {
      val VTP = variantTypeIsPoLOF(vtype);
      val NS  = txIsLofSensitive(tx)
      if(NS && VTP){
        (true,Seq(tx+":"+vtype+":"+aa+":LS"));
      } else if(VTP){
        (false,Seq(tx+":"+vtype+":"+aa+":NLS"));
      } else {
        (false,Seq[String]());
      }
    }}.foldLeft( ACMGCrit("PVS",1,false,Seq[String]()) ){case (soFar,(isFlagged,opStr)) => {
      ACMGCrit("PVS",1,isFlagged,opStr).merge(soFar);
    }}
  }
  
  def variantTypeIsPoLOF(vtype : Array[String]) : Boolean = {
    vtype(1) == "LLOF" || vtype(1) == "PLOF";
  }
  
  
  def getAttributeDoubleListOption(v : VariantContext, tag : String, ct : Int, defaultVal : Double = -1.0) : Vector[Double] = {
    if(v.hasAttribute(tag)) v.getAttributeAsList(tag).toVector.map(x => if(x.toString == ".") defaultVal 
                                                                        else string2double(x.toString)) 
    else repToSeq(defaultVal,ct).toVector;
  }
  def checkDatabaseAlleleFreq(af : Vector[Double], altIdx : Int, tag : String, soFar : ACMGCrit, thresh : Double = 0.05) : (Boolean,ACMGCrit) = {
    //val attrAF = getAttributeDoubleListOption(v,tag,0);
        if(af(altIdx) > thresh){
          return (true, ACMGCrit("BA",1,true,Vector(tag+"="+af(altIdx))).merge(soFar))
        } else {
          return (true, ACMGCrit("BA",1,false,Vector(tag+"="+af(altIdx))).merge(soFar))
        }
  }
      
  def checkDbAlleFreq(v : VariantContext, altIdx : Int, thresh : Double = 0.05) : (ACMGCrit,ACMGCrit) = {
    val TKG = getAttributeDoubleListOption(v,"1KG_AF",0);
    val ESP = getAttributeDoubleListOption(v,"ESP_EA_AF",0);
    val EXAC = getAttributeDoubleListOption(v,"dbNSFP_ExAC_AF",0);
    
    //var out : Option[ACMGCrit] = None;
    
    val (tkgFlag,step0)  = checkDatabaseAlleleFreq(TKG,altIdx,"1KG_AF",ACMGCrit("BA",1,false,Seq[String]()),thresh=thresh);
    val (espFlag,step1)  = checkDatabaseAlleleFreq(ESP,altIdx,"ESP_EA_AF",step0,thresh=thresh);
    val (exacFlag,step2) = checkDatabaseAlleleFreq(EXAC,altIdx,"dbNSFP_ExAC_AF",step1,thresh=thresh);
    
    val PM1 = if(TKG(altIdx) <= 0 && ESP(altIdx) <= 0 && EXAC(altIdx) <= 0){
      ACMGCrit("PM",1,true,Vector())
    } else {
      ACMGCrit("PM",1,false,Vector())
    }
    return (PM1,step2);
  }
  
  case class ProteinVariant(tx : String, varType : String, start : Int, end : Int, ref : String, alt : String){
    
  }
  /*case class ProteinVariantDEL(tx : String, pos : Int, end:Int) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = false;
  }
  case class ProteinVariantMIS(tx : String, pos : Int, alt : String) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = false;
  }
  case class ProteinVariantINS(tx : String, pos : Int, alt : String) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = false;
  }
  case class ProteinVariantFS(tx : String, pos : Int, alt :String) extends ProteinVariant(tx,pos){
    def isSevere : Boolean = true;
  }
  case class ProteinVariantDELINS(tx : String, pos : Int, end : Int, alt : String) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = false;
  }
  case class ProteinVariantSPLICE(tx : String, pos : Int, end : Int, alt : String) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = true;
  }
  case class ProteinVariantSTOPLOSS(tx : String, pos : Int) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = false;
  }
  case class ProteinVariantSTARTLOSS(tx : String, pos : Int) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = true;
  }
  case class ProteinVariantNONS(tx : String, pos : Int) extends ProteinVariant(tx,pos) {
    def isSevere : Boolean = true;
  }*/
  
  
  /*
   * Bottom-level output commands:
   */
  
  case class ACMGCrit(x : String, y : Int, flag : Boolean, anno : Seq[String]) {
    def checkValid : Boolean = {
      return true;
    }
    def merge(c : ACMGCrit) : ACMGCrit = {
      if(c.x != x && c.y != y) error("Error! Attempt to merge unequal ACMGCrit!");
      return ACMGCrit(x,y,flag || c.flag, c.anno ++ this.anno);
    }
  }
  
  implicit object ACMGCritOrdering extends Ordering[ACMGCrit] {
    val ord = Vector("BP","BS","BA","PP","PM","PS","PVS").reverse;
    def compare(x : ACMGCrit, y : ACMGCrit) : Int = {
      val xx = ord.indexOf(x.x);
      val yx = ord.indexOf(y.x);
      if(xx == yx){
        return implicitly[Ordering[Int]].compare(xx,yx);
      } else {
        return implicitly[Ordering[Int]].compare(x.y,y.y);
      }
    }
  }
  
  class ACMGCritSet() {
    var critSet : Set[ACMGCrit] = Set[ACMGCrit]();
    def addCrit(c : ACMGCrit){
      critSet = critSet + c;
    }
    def addCrit(c : Option[ACMGCrit]){
      c match {
        case Some(c) => {
          addCrit(c);
        }
        case None => {
          //do nothing.
        }
      }
    }
    def getCritSet : Seq[ACMGCrit] = critSet.toVector.sorted;
  }


  
  def getACMGPathogenicityRating(inputCriteria : Seq[ACMGCrit]) : String = {
    val criteria = inputCriteria.filter(_.flag);
    val numPVS = criteria.count{(c) => c.x == "PVS"}
    val numPS  = criteria.count{(c) => c.x == "PS"}
    val numPM  = criteria.count{(c) => c.x == "PM"}
    val numPP  = criteria.count{(c) => c.x == "PP"}
    val numBA  = criteria.count{(c) => c.x == "BA"}
    val numBS  = criteria.count{(c) => c.x == "BS"}
    val numBP  = criteria.count{(c) => c.x == "BP"}
    
    val isPathogenic = (numPVS >= 1 && (
                           (numPS >= 1) || 
                           (numPM >= 2) ||
                           (numPM + numPP >= 2) ||
                           (numPP >= 2)
                       )) || 
                       ( numPS >= 2 ) ||
                       ( numPS >= 1 && (
                           (numPM >= 3) ||
                           (numPM >= 2 && numPP >= 2) ||
                           (numPM >= 1 && numPP >= 4)
                       ));
    
    val isLikelyPatho = (numPVS >= 1 && numPM >= 1) ||
                        (numPS  >= 1 && numPM >= 1) ||
                        (numPS  >= 1 && numPP >= 2) || //CHECK ME???
                        (numPM  >= 3) ||
                        (numPM  >= 2 && numPP >= 2) ||
                        (numPM  >= 1 && numPP >= 4);
    
    val isLikelyBenign = (numBS >= 1 && numBP >= 1) ||
                         (numBP >= 2);
    val isBenign  = (numBA >= 1) ||
                    (numBS >= 2);
    
    if(isPathogenic){
      if(isBenign){
        return "VUS";
      } else {
        return "PATHO";
      }
    }
    if(isBenign){
      return("B");
    }
    if(isLikelyPatho){
      if(isLikelyBenign){
        return("VUS");
      } else {
        return("LPATH");
      }
    }
    if(isLikelyBenign){
      return("LB");
    }
    return("VUS")
  }
  
}









