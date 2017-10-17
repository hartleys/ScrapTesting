package internalUtils

import htsjdk.variant._;
import htsjdk.variant.variantcontext._;
import htsjdk.variant.vcf._;
import java.io.File;
import scala.collection.JavaConversions._
import java.io._;
import internalUtils.commandLineUI._;
import internalUtils.Reporter._;

import scala.collection.JavaConversions._
import scala.collection.JavaConverters._

import internalUtils.optionHolder._;
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import internalUtils.commandLineUI._;
import internalUtils.fileUtils._;
import internalUtils.TXUtil._;

import htsjdk.variant._;
import htsjdk.variant.variantcontext._;
import htsjdk.variant.vcf._;

import internalUtils.genomicUtils._;
import internalUtils.commonSeqUtils._;


object VcfTool {
  
  val TOP_LEVEL_VCF_TAG : String = "SWH_";
  
  case class VCFAnnoCodes(
        txList_TAG : String = TOP_LEVEL_VCF_TAG+ "TXLIST",
        vType_TAG : String = TOP_LEVEL_VCF_TAG+"VARTYPE",
        vMutG_TAG : String = TOP_LEVEL_VCF_TAG+"varG",
        vMutR_TAG : String = TOP_LEVEL_VCF_TAG+"varR",
        vMutC_TAG : String = TOP_LEVEL_VCF_TAG+"varC",
        vMutP_TAG : String = TOP_LEVEL_VCF_TAG+"varPredP",
        
        vTypeShort_TAG : String = TOP_LEVEL_VCF_TAG+"ShortVARTYPE",
        vMutPShort_TAG : String = TOP_LEVEL_VCF_TAG+"ShortVarPredP",
        vMutLVL_TAG : String = TOP_LEVEL_VCF_TAG + "varLVL",
        vMutINFO_TAG : String = TOP_LEVEL_VCF_TAG + "varRAWDATA",
        
        grpAC_TAG : String = TOP_LEVEL_VCF_TAG+"AC_GRP_",
        grpAF_TAG : String = TOP_LEVEL_VCF_TAG+"AF_GRP_",
        grpHomCt_TAG : String = TOP_LEVEL_VCF_TAG+"HomCt_GRP_",
        grpHetCt_TAG : String = TOP_LEVEL_VCF_TAG+"HetCt_GRP_",
        grpHomFrq_TAG : String = TOP_LEVEL_VCF_TAG+"HomFrq_GRP_",
        grpHetFrq_TAG : String = TOP_LEVEL_VCF_TAG+"HetFrq_GRP_",
        grpMisCt_TAG : String = TOP_LEVEL_VCF_TAG+"MisCt_GRP_",
        grpMisFrq_TAG : String = TOP_LEVEL_VCF_TAG+"MisFrq_GRP_",
        grpRefAC_TAG : String = TOP_LEVEL_VCF_TAG+"RefAC_GRP_",
        grpRefAF_TAG : String = TOP_LEVEL_VCF_TAG+"RefAF_GRP_",
        
        isSplitMulti_TAG : String = TOP_LEVEL_VCF_TAG+"isSplitMult",
        splitIdx_TAG     : String = TOP_LEVEL_VCF_TAG+"splitIdx",
        numSplit_TAG     : String = TOP_LEVEL_VCF_TAG+"numSplit",
        splitAlle_TAG     : String = TOP_LEVEL_VCF_TAG+"fullAlleList",
        
        assess_IsRepetitive : String = TOP_LEVEL_VCF_TAG+"locusIsRep",
        assess_IsConserved : String = TOP_LEVEL_VCF_TAG+"locusIsCons",
        assess_IsHotspot : String = TOP_LEVEL_VCF_TAG+"locusIsHotspot",
        assess_domain : String = TOP_LEVEL_VCF_TAG+"locusDomain",
        
        assess_criteria : String = TOP_LEVEL_VCF_TAG+"ACMG_criteria",
        
        assess_GeneSenseLOF :String = TOP_LEVEL_VCF_TAG+"ACMG_GeneSenseLOF",
        assess_GeneSenseMis :String = TOP_LEVEL_VCF_TAG+"ACMG_GeneSenseMis",
        assess_ctrlAFMIN : String = TOP_LEVEL_VCF_TAG+"ACMG_ctrlAFMIN",
        assess_ctrlAFMAX : String = TOP_LEVEL_VCF_TAG+"ACMG_ctrlAFMAX",
        
        assess_LofGenes : String =  TOP_LEVEL_VCF_TAG+"ACMG_LofGenes",
        assess_MisGenes : String =  TOP_LEVEL_VCF_TAG+"ACMG_MisGenes",
        
        assess_PVS1 : String = TOP_LEVEL_VCF_TAG+"ACMG_PVS1",
        assess_PS1 : String = TOP_LEVEL_VCF_TAG+"ACMG_PS1",
        assess_PM1 : String = TOP_LEVEL_VCF_TAG+"ACMG_PM1",
        assess_PM2 : String = TOP_LEVEL_VCF_TAG+"ACMG_PM2",
        assess_PM4 : String = TOP_LEVEL_VCF_TAG+"ACMG_PM4",
        assess_PM5 : String = TOP_LEVEL_VCF_TAG+"ACMG_PM5",
        assess_PP2 : String = TOP_LEVEL_VCF_TAG+"ACMG_PP2",
        assess_PP3 : String = TOP_LEVEL_VCF_TAG+"ACMG_PP3",
        assess_BP1 : String = TOP_LEVEL_VCF_TAG+"ACMG_BP1",
        assess_BP2 : String = TOP_LEVEL_VCF_TAG+"ACMG_BP2",
        assess_BP3 : String = TOP_LEVEL_VCF_TAG+"ACMG_BP3",
        assess_BP4 : String = TOP_LEVEL_VCF_TAG+"ACMG_BP4",
        assess_BP7 : String = TOP_LEVEL_VCF_TAG+"ACMG_BP7",
        assess_BS1 : String = TOP_LEVEL_VCF_TAG+"ACMG_BS1",
        assess_BS2 : String = TOP_LEVEL_VCF_TAG+"ACMG_BS2",
        assess_BA1 : String = TOP_LEVEL_VCF_TAG+"ACMG_BA1",
        assess_RATING : String = TOP_LEVEL_VCF_TAG+"ACMG_RATING",
        assess_WARNINGS : String = TOP_LEVEL_VCF_TAG+"ACMG_WARN",
        assess_WARNFLAG : String = TOP_LEVEL_VCF_TAG+"ACMG_WARNFLAG",
        
        assess_exactMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_ExactMatchRS",
        assess_aminoMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_AminoMatchRS",
        assess_nearMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_NearMatchRS",
        
        geneIDs : String = TOP_LEVEL_VCF_TAG+"txGeneIDs",
        
        delims : Vector[String] = Vector(",","|")
      );
  
  val UNAUTOMATED_ACMG_PARAMS = Seq[(String,Int,String)](
        ("PS",2 , TOP_LEVEL_VCF_TAG+"ACMG_PS2"),
        ("PS",3 , TOP_LEVEL_VCF_TAG+"ACMG_PS3"),
        ("PS",4 , TOP_LEVEL_VCF_TAG+"ACMG_PS4"),
        ("PM",3 , TOP_LEVEL_VCF_TAG+"ACMG_PM3"),
        ("PM",6, TOP_LEVEL_VCF_TAG+"ACMG_PM6"),
        ("PP",1 , TOP_LEVEL_VCF_TAG+"ACMG_PP1"),
        ("PP",4 , TOP_LEVEL_VCF_TAG+"ACMG_PP4"),
        ("PP",5 , TOP_LEVEL_VCF_TAG+"ACMG_PP5"),
        ("BP",2 , TOP_LEVEL_VCF_TAG+"ACMG_BP2"),
        ("BP",5 , TOP_LEVEL_VCF_TAG+"ACMG_BP5"),
        ("BP",6 , TOP_LEVEL_VCF_TAG+"ACMG_BP6"),
        ("BS",1 , TOP_LEVEL_VCF_TAG+"ACMG_BS1"),
        ("BS",2 , TOP_LEVEL_VCF_TAG+"ACMG_BS2"),
        ("BS",3 , TOP_LEVEL_VCF_TAG+"ACMG_BS3"),
        ("BS",4 , TOP_LEVEL_VCF_TAG+"ACMG_BS4")
      )
  
    //  val isRepetitive = locusIsRepetitive(v.getContig(),v.getStart());
    //val isConserved  = locusIsConserved(v.getContig(),v.getStart());
    //val isHotspot    = locusIsHotspot(v.getContig(),v.getStart());
    
  abstract class VCFWalker{
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader);
    
    def walkVCFFile(infile : String, outfile : String, chromList : Option[List[String]], vcfCodes : VCFAnnoCodes = VCFAnnoCodes(), verbose : Boolean = true){
      val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(infile, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes);
      val (vcIter2, newHeader) = this.walkVCF(vcIter = vcIter, vcfHeader = vcfHeader)
    
      val vcfWriter = internalUtils.VcfTool.getVcfWriter(outfile, header = newHeader);
    
      vcIter2.foreach(vc => {
        vcfWriter.add(vc)
      })
      vcfWriter.close();
    }
    
  }
  
    def getAltAllelesInOrder(vc : VariantContext) : Seq[Allele] = {
      Range(0,vc.getNAlleles()-1).map(i => {
        vc.getAlternateAllele(i);
      })
    }
    def getAllelesInOrder(vc : VariantContext) : Seq[Allele] = {
      vc.getReference() +: getAltAllelesInOrder(vc);
    }
  
  //htsjdk.variant.variantcontext.writer.VariantContextWriter
  
  def addHeaderLines(vcfHeader : VCFHeader, headerLines : Seq[VCFHeaderLine]) : VCFHeader = {
    val newHeaderLines : List[VCFHeaderLine] = vcfHeader.getMetaDataInInputOrder().toList ++
                                               headerLines.toList;
    
    
    
    return new VCFHeader(newHeaderLines.toSet.asJava,
                         vcfHeader.getSampleNamesInOrder().asInstanceOf[java.util.List[String]]);
  }
  def replaceHeaderLines(vcfHeader : VCFHeader, headerLines : Seq[VCFHeaderLine]) : VCFHeader = {
    val newHeaderLines : List[VCFHeaderLine] = headerLines.toList;

    return new VCFHeader(newHeaderLines.toSet.asJava,
                         vcfHeader.getSampleNamesInOrder().asInstanceOf[java.util.List[String]]);
  }
  
  
  def getVcfWriter(outfile : String, header : VCFHeader) : htsjdk.variant.variantcontext.writer.VariantContextWriter = {
    val vcfb = new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder();
    
    //vcfb.setReferenceDictionary(vcfHeader.getSequenceDictionary());
    vcfb.unsetOption(htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY);
    vcfb.setOutputFile(outfile);
    
    if(outfile.takeRight(3) == ".gz"){
      vcfb.setOutputFileType(htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
    } else {
      vcfb.setOutputFileType(htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF);
    }
        
    val vcfWriter : htsjdk.variant.variantcontext.writer.VariantContextWriter = vcfb.build();
    
    vcfWriter.writeHeader(header);
    
    return vcfWriter;
  }
  
  def getVcfIterator(   infile : String, 
                        chromList : Option[List[String]],
                        vcfCodes : VCFAnnoCodes = VCFAnnoCodes(),
                        progressVerbosity : (Int,Int,Int) = (10000,50000,100000)
                     ) : (Iterator[VariantContext], VCFHeader) = {
      val chromSet = chromList match {
        case Some(lst) => Some(lst.toSet);
        case None => None;
      }
      reportln("Starting VCF read...","progress");
      
      val vcfReader = new VCFFileReader(new File(infile),false);
      val vcfHeader = vcfReader.getFileHeader();
      val vcfIterator : Iterator[VariantContext] = vcfReader.iterator(); //chromSet match {
      val (a,b,c) = progressVerbosity;
      
      val finalIterator = chromSet match {
        case Some(cs) => {
          internalUtils.stdUtils.wrapIteratorWithAdvancedProgressReporter[VariantContext](vcfIterator,
              internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevel[VariantContext](elementTitle = "lines", a,b,c,((a : VariantContext,i : Int) => a.getContig()))).filter(v =>  cs.contains(v.getContig()));
        }
        case None => {
          internalUtils.stdUtils.wrapIteratorWithAdvancedProgressReporter[VariantContext](vcfIterator,
              internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[VariantContext](elementTitle = "lines", lineSec = 60))
        }
      }
      
      //val wrappedVcfIterator = progressVerbosity match {
      //  //case Some((a,b,c)) => internalUtils.stdUtils.wrapIteratorWithProgressReporter(vcfIterator , internalUtils.stdUtils.IteratorProgressReporter_ThreeLevel("lines",a,b,c) )
      //  case Some((a,b,c)) => internalUtils.stdUtils.wrapIteratorWithAdvancedProgressReporter[VariantContext](filteredIterator,
      //                           internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[VariantContext](elementTitle = "lines", lineSec = 60))
      //                        //internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAccel[VariantContext](elementTitle = "lines", 
      //                        //                                   accelFactor  = 10,
      //                        //                                   maxAccel = a,
      //                        //                                   dotThreshold = 1,
      //                        //                                   dotSpaceThreshold = 5,
      //                        //                                  dotNewlineThreshold = 10))
      //  case None => vcfIterator
      //}
      //val finalIterator = chromSet match {
      //      case Some(chroms) => wrappedVcfIterator.filter(v =>  chroms.contains(v.getContig()));
      //      case None => wrappedVcfIterator;
      //}
      return (finalIterator,vcfHeader);
  }

  
  
/*
   abstract class VcfMetaLine {
     def key : String;
     def value : String;
     override def toString() : String = "##"+key+"="+value;
   }
   case class VcfAnnoLine(k : String, v : String) extends VcfMetaLine {
     def key : String = k;
     def value : String = v;
   }
   case class VcfInfoLine(ID : String, Number : String, Type : String, Description : String, Source:Option[String] = None, Version:Option[String] = None) extends VcfMetaLine {
     def key : String = "INFO";
     def value : String = "<" + 
                          "ID="+ ID+""+
                          ",Number="+Number+""+
                          ",Type="+Type+""+
                          ",Description="+Description+
                          (if(Source.isEmpty) "" else ",Source="+Source.get)+
                          (if(Version.isEmpty) "" else ",Version="+Version.get)+
                          ">";
   }
   case class VcfFilterLine(ID : String, Description : String) extends VcfMetaLine {
     def key : String = "FILTER";
     def value : String = "<"+
                          "ID="+ ID+""+
                          ",Description="+Description+
                          ">";
   }
   case class VcfFormatLine(ID : String, Number : String, Type : String, Description : String) extends VcfMetaLine {
     def key : String = "FORMAT";
     def value : String = "<" + 
                          "ID="+ ID+""+
                          ",Number="+Number+""+
                          ",Type="+Type+""+
                          ",Description="+Description+
                          ">";
   }
   
   private def parseVcfMetadataValue(s : String) : Seq[(String,String)] = {
     if(s.length < 2 || s.head != '<' || s.last != '>'){
       error("FATAL ERROR: Malformed metadata line in VCF file: value is too short or isn't bound by angle-brackets (errCode VcfTool:parseVcfMetadataValue:65)!");
     }
     val valCells = internalUtils.stdUtils.parseTokens(s.init.tail,',');
     val valPairs = valCells.map((v) => {
       val pair = internalUtils.stdUtils.parseTokens(v,'=');
       if(pair.length != 2) error("FATAL ERROR: Malformed metadata line in VCF file (errCode VcfTool:76)!");
       (pair(0),pair(1));
     });
     return valPairs;
   }
   
   def readVcfMetadata(lines : Iterable[String], sampleLine : String) : VcfMetadata = {
     val metaLines : Seq[VcfMetaLine] = lines.map((s : String) => readVcfMetaLine(s)).toSeq;
     val (info : Seq[VcfInfoLine], filter : Seq[VcfFilterLine], format : Seq[VcfFormatLine], anno : Seq[VcfAnnoLine]) = {
       metaLines.foldLeft((Seq[VcfInfoLine](),Seq[VcfFilterLine](),Seq[VcfFormatLine](),Seq[VcfAnnoLine]()))((soFar,curr) => {
         curr match {
           case x : VcfInfoLine => (soFar._1 :+ x, soFar._2, soFar._3, soFar._4);
           case x : VcfFilterLine => (soFar._1, soFar._2 :+ x, soFar._3, soFar._4);
           case x : VcfFormatLine => (soFar._1 , soFar._2, soFar._3 :+ x, soFar._4);
           case x : VcfAnnoLine => (soFar._1, soFar._2, soFar._3, soFar._4 :+ x);
         }
       }) 
     }
     val sampleCells = sampleLine.split("\t");
     if(sampleCells.length < 9) error("FATAL ERROR: Malformed Vcf Header line. Less than 9 columns! (errCode VcfTool:readVcfMetadata:72)\n offending line:\""+sampleLine+"\"");
     val sampleID = sampleCells.slice(9,sampleCells.length);
     
     return VcfMetadata(info, filter, format, anno, sampleID);
   }
   
   def readVcfMetaLine(line : String) : VcfMetaLine = {
     if(line.substring(0,2) != "##"){
       error("FATAL ERROR: Impossible state! readVcfMetaLine has been given a line that does not start with \"##\"! (errcode VcfTool:readVcfMetaLine:78)");
     }
     val cells = line.substring(2).split("=",2);
     if(cells(0) == "INFO"){
       readVcfInfoLine(cells(1));
     } else if(cells(0) == "FILTER"){
       readVcfFilterLine(cells(1));
     } else if(cells(0) == "FORMAT"){
       readVcfFormatLine(cells(1));
     } else {
       VcfAnnoLine(cells(0), cells(1));
     }
   }
   private def readVcfInfoLine(v : String) : VcfInfoLine = {
       val valMap = parseVcfMetadataValue(v).toMap;
       if(! valMap.contains("ID")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no ID key! (errCode VcfTool:readVcfInfoLine:84)");
       if(! valMap.contains("Number")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no Number key! (errCode VcfTool:readVcfInfoLine:85)");
       if(! valMap.contains("Type")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no Type key! (errCode VcfTool:readVcfInfoLine:86)");
       if(! valMap.contains("Description")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no Description key! (errCode VcfTool:readVcfInfoLine:87)");
       return VcfInfoLine(valMap("ID"), valMap("Number"), valMap("Type"), valMap("Description"), valMap.get("Source"), valMap.get("Version"));
   }
   private def readVcfFilterLine(v : String) : VcfFilterLine = {
       val valMap = parseVcfMetadataValue(v).toMap;
       if(! valMap.contains("ID")) error("FATAL ERROR: Malformed FILTER metadata line in VCF file: no ID key! (errCode VcfTool:readVcfFilterLine:84)");
       if(! valMap.contains("Description")) error("FATAL ERROR: Malformed FILTER metadata line in VCF file: no Description key! (errCode VcfTool:readVcfFilterLine:87)");
       return VcfFilterLine(valMap("ID"), valMap("Description"));
   }
   private def readVcfFormatLine(v : String) : VcfFormatLine = {
       val valMap = parseVcfMetadataValue(v).toMap;
       if(! valMap.contains("ID")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no ID key! (errCode VcfTool:readVcfFormatLine:84)");
       if(! valMap.contains("Number")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no Number key! (errCode VcfTool:readVcfFormatLine:85)");
       if(! valMap.contains("Type")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no Type key! (errCode VcfTool:readVcfFormatLine:86)");
       if(! valMap.contains("Description")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no Description key! (errCode VcfTool:readVcfFormatLine:87)");
       return VcfFormatLine(valMap("ID"), valMap("Number"), valMap("Type"), valMap("Description"));
   }
   
   case class VcfMetadata(info : Seq[VcfInfoLine],filter : Seq[VcfFilterLine], format : Seq[VcfFormatLine], anno : Seq[VcfAnnoLine], sampleID : Seq[String]) {
     lazy val metaLines   : Seq[VcfMetaLine] = info ++ filter ++ format ++ anno;

     lazy val infoMap : Map[String,VcfInfoLine] = info.map((v) =>{
       (v.ID,v);
     }).toMap;
     lazy val filterMap : Map[String,VcfFilterLine] = filter.map((v) =>{
       (v.ID,v);
     }).toMap;
     lazy val formatMap : Map[String,VcfFormatLine] = format.map((v) =>{
       (v.ID,v);
     }).toMap;
     
   }
   
   abstract class VcfLine {
     def CHROM : String;
     def POS : Int;
     def ID : String;
     def REF : String;
     def ALT: String;
     def QUAL : Double;
     def FILTER: String;
     def INFO: String;
     def FORMAT : String;
     def GENOTYPES: Seq[String];
     def metadata : VcfMetadata;
     
     override def toString() : String = CHROM +"\t"+POS+"\t"+ID+"\t"+REF+"\t"+ALT+"\t"+QUAL+"\t"+FILTER+"\t"+INFO+"\t"+FORMAT+"\t"+GENOTYPES.mkString("\t");
     
     def fmt : Seq[String];

     lazy val idSeq : Seq[String] = ID.split(";");
     lazy val altSeq : Seq[String] = ALT.split(",");
     
     lazy val fmtMeta : Seq[VcfFormatLine] = fmt.map((f : String) =>{
       metadata.formatMap.get(f) match {
         case Some(x) => x;
         case None => {
           error("FATAL VCF PARSING ERROR: FORMAT ID not found!\n offending line:\""+toString()+"\"");
           null;
         }
       }
     })
     lazy val fmtInfo : Seq[(String,String)] = fmtMeta.map((m : VcfFormatLine)=> (m.Number,m.Type));
   }
   
   case class InputVcfLine(line : String, meta : VcfMetadata) extends VcfLine {
     def metadata = meta;
     lazy val cells = {
       val c = line.split("\t");
       if(c.length < 9) error("FATAL ERROR: Vcf Line has fewer than 9 columns:\n  Offending line: \""+line+"\"");
       c;
     }
     def CHROM : String = cells(0);
     lazy val pos : Int = internalUtils.stdUtils.string2int(cells(1));
     def POS : Int = pos;
     def ID : String = cells(2);
     def REF : String = cells(3);
     def ALT : String = cells(4);
     lazy val qual : Double = internalUtils.stdUtils.string2double(cells(5));
     def QUAL : Double = qual;
     def FILTER : String = cells(6);
     def INFO : String = cells(7);
     def FORMAT : String = cells(8);
     def GENOTYPES : Seq[String] = cells.slice(9,cells.length);
     
     //def passFilter : Boolean = 
     
     lazy val FMT : Seq[String] = FORMAT.split(":");
     def fmt : Seq[String] = FMT;
     
     lazy val genoTableBySample : Seq[Seq[String]] = GENOTYPES.map((g : String) => {
       val out = g.split(":").toSeq;
       if(out.length != fmt.length) error("FATAL ERROR: Vcf genotype has the wrong number of elements.\n offending line: \""+line+"\"");
       out;
     });
     lazy val genoTable : Seq[Seq[String]] = internalUtils.stdUtils.transposeMatrix(genoTableBySample);
   }
   //case class OutputVcfLine(chrom : String, pos : Int, id : String, ref : String, alt : String, qual : Double, filter : String, Info : String, genotypes : Seq[String]){
     
   //}
   
   /*
   abstract class VcfMetadata {
     def getMetaLines   : Seq[VcfMetaLine];
     def getInfoLines   : Seq[VcfInfoLine];
     def getFilterLines : Seq[VcfFilterLine];
     def getFormatLines : Seq[VcfFormatLine];
     def getInfoMap : Map[String,VcfInfoLine] = getInfoLines.map((v) =>{
       (v.ID,v);
     }).toMap;
     def getFormatMap : Map[String,VcfFormatLine] = getFormatLines.map((v) =>{
       (v.ID,v);
     }).toMap;
     def getFilterMap : Map[String,VcfFilterLine] = getFilterLines.map((v) =>{
       (v.ID,v);
     }).toMap;
   }
   */
   
   abstract class VcfNumberField {
     def isInteger : Boolean;
     def isSpecial : Boolean;
     def isKnown : Boolean;
   }
   case class VcfNumberFieldDot() extends VcfNumberField {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def isKnown : Boolean = false;
   }
   abstract class VcfNumberFieldKnown extends VcfNumberField {
     def isKnown : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int;
   }
   case class VcfNumberFieldInt(value : Int) extends VcfNumberFieldKnown {
     def isInteger : Boolean = true;
     def isSpecial : Boolean = false;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = value;
   }
   case class VcfNumberFieldA() extends VcfNumberFieldKnown {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = altAlleleCt;
   }
   case class VcfNumberFieldR() extends VcfNumberFieldKnown {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = altAlleleCt + 1;
   }
   case class VcfNumberFieldG() extends VcfNumberFieldKnown {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = genotypeCt;
   }
   
   abstract class VcfVal {
     def isInt : Boolean = false;
     def isFloat : Boolean = false;
     def isChar : Boolean = false;
     def isString : Boolean = false;
   }
   case class VcfInt(v : Int) extends VcfVal {
     override def toString() = v.toString();
     override def isInt : Boolean = true;
   }
   case class VcfFloat(v : Double) extends VcfVal {
     override def toString() = v.toString();
     override def isFloat : Boolean = true;
   }
   case class VcfChar(v : Char) extends VcfVal {
     override def toString() = v.toString();
     override def isChar : Boolean = true;
   }
   case class VcfString(v : String) extends VcfVal {
     override def toString() = v.toString();
     override def isString : Boolean = true;
   }
   
   def readVcfVal(v : String, fmt : String) : VcfVal = {
     if(fmt == "Integer"){
       return VcfInt(internalUtils.stdUtils.string2int(v));
     } else if(fmt == "Float"){
       return VcfFloat(internalUtils.stdUtils.string2double(v));
     } else if(fmt == "Char"){
       if(v.length > 1) error("FATAL ERROR: Malformed VCF. Found 'char' formatted field with length > 1!");
       return VcfChar(v.charAt(0));
     } else if(fmt == "String"){
       return VcfString(v);
     } else {
       error("FATAL ERROR: Malformed VCF. Unrecognized format: \""+fmt+"\"");
       return null;
     }
   } 
   
   
   def getVcfReader(lines : Iterator[String]) : (VcfMetadata, Iterator[InputVcfLine]) = {
     val (metaLines, bodyLines) = internalUtils.stdUtils.splitIterator(lines, (ln : String) => {
       ln.startsWith("##");
     });
     if(! bodyLines.hasNext) error("FATAL ERROR: Vcf File does not have header or body lines.");
     val headerLine = bodyLines.next;
     if(! headerLine.startsWith("#")) error("FATAL ERROR: Vcf File header line not present or does not start with \"#\"");
     if(! bodyLines.hasNext) error("FATAL ERROR: Vcf File does not have body lines.");
     
     val meta : VcfMetadata = readVcfMetadata(metaLines,headerLine);
     
     val vcfLines : Iterator[InputVcfLine] = bodyLines.map( (line : String) => {
       InputVcfLine(line,meta);
     });
     
     return ((meta,vcfLines));
   }
   
  */
}















