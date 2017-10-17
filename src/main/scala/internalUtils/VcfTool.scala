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
  
  def getAttributeAsStringList(vc : VariantContext, tag : String) : Seq[String] = {
    val lst = vc.getAttributeAsList(tag).asScala;
    if(lst.length == 1){
      return( lst.head.toString.split(",") );
    } else {
      return( lst.map(_.toString()) );
    }
  }
  
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
        grpRefAC_TAG : String = TOP_LEVEL_VCF_TAG+"AC_REF_GRP_",
        grpRefAF_TAG : String = TOP_LEVEL_VCF_TAG+"AF_REF_GRP_",
        
        grpAltAC_TAG : String = TOP_LEVEL_VCF_TAG+"AC_ALT_GRP_",
        grpAltHomCt_TAG : String = TOP_LEVEL_VCF_TAG+"HomCt_ALT_GRP_",
        grpAltHetCt_TAG : String = TOP_LEVEL_VCF_TAG+"HetCt_ALT_GRP_",
        grpRefHomCt_TAG : String = TOP_LEVEL_VCF_TAG+"HomCt_REF_GRP_",
        grpRefHetCt_TAG : String = TOP_LEVEL_VCF_TAG+"HetCt_REF_GRP_",
        
        CLNVAR_SUMSIG : String = TOP_LEVEL_VCF_TAG+"ClinVar_SumSig",
        CLNVAR_SUMSIGWARN : String = TOP_LEVEL_VCF_TAG+"ClinVar_SumSig_Warning",
        CLNVAR_SUMSIG_SAFE : String = TOP_LEVEL_VCF_TAG+"ClinVar_SumSig_NoWarn",
        
        domainIds : String = TOP_LEVEL_VCF_TAG+"domainIdList",
        
        isSplitMulti_TAG : String = TOP_LEVEL_VCF_TAG+"isSplitMult",
        splitIdx_TAG     : String = TOP_LEVEL_VCF_TAG+"splitIdx",
        numSplit_TAG     : String = TOP_LEVEL_VCF_TAG+"numSplit",
        splitAlle_TAG     : String = TOP_LEVEL_VCF_TAG+"fullAlleList",
        
        assess_IsRepetitive : String = TOP_LEVEL_VCF_TAG+"locusIsRep",
        assess_IsConserved : String = TOP_LEVEL_VCF_TAG+"locusIsCons",
        assess_IsHotspot : String = TOP_LEVEL_VCF_TAG+"locusIsHotspot",
        assess_domain : String = TOP_LEVEL_VCF_TAG+"locusDomain",
        assess_IsMappable : String = TOP_LEVEL_VCF_TAG+"locusIsMappable",
        assess_pseudogene : String = TOP_LEVEL_VCF_TAG+"locusIsPseudogene", 

        assess_criteria : String = TOP_LEVEL_VCF_TAG+"ACMG_criteria",
        
        assess_GeneSenseLOF :String = TOP_LEVEL_VCF_TAG+"ACMG_GeneSenseLOF",
        assess_GeneSenseMis :String = TOP_LEVEL_VCF_TAG+"ACMG_GeneSenseMis",
        assess_ctrlAFMIN : String = TOP_LEVEL_VCF_TAG+"ACMG_ctrlAFMIN",
        assess_ctrlAFMAX : String = TOP_LEVEL_VCF_TAG+"ACMG_ctrlAFMAX",
        
        assess_geneList : String =  TOP_LEVEL_VCF_TAG+"ACMG_genes",
        assess_geneTxList : String =  TOP_LEVEL_VCF_TAG+"ACMG_geneTX",
        
        assess_LofGenes : String =  TOP_LEVEL_VCF_TAG+"ACMG_LofGenes",
        assess_MisGenes : String =  TOP_LEVEL_VCF_TAG+"ACMG_MisGenes",
        
        assess_LofTX : String =  TOP_LEVEL_VCF_TAG+"ACMG_LofTX",
        assess_MisTX : String =  TOP_LEVEL_VCF_TAG+"ACMG_MisTX",
        
        assess_geneLofTxRatio : String =  TOP_LEVEL_VCF_TAG+"ACMG_geneLofTxRatio",
        assess_geneMisTxRatio : String =  TOP_LEVEL_VCF_TAG+"ACMG_geneMisTxRatio",
        
        assess_refSeqLOF : String = TOP_LEVEL_VCF_TAG+"ACMG_geneCanonLOF",
        assess_refSeqMis : String = TOP_LEVEL_VCF_TAG+"ACMG_geneCanonMis",
        assess_refSeqKnown : String = TOP_LEVEL_VCF_TAG+"ACMG_geneCanonIsKnown",
        
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
        
        assess_PVS1_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_PVS1_CANON",
        assess_PP2_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_PP2_CANON",
        assess_BP1_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_BP1_CANON",
        assess_BP3_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_BP3_CANON",
        assess_PM4_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_PM4_CANON",
        assess_PS1_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_PS1_CANON",
        assess_PM5_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_PM5_CANON",
        assess_BP7_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_BP7_CANON",
        assess_RATING_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_RATING_CANON", 
        
        //NEW: 
         
        assess_forBT_codonLOF : String = TOP_LEVEL_VCF_TAG+"ASSESS_LOF_codonPathoRS",
        assess_forBT_downstreamLOF  : String = TOP_LEVEL_VCF_TAG+"ASSESS_LOF_downstreamPathoRS",
        
        assess_forBT_codonLOF_CANON : String = TOP_LEVEL_VCF_TAG+"ASSESS_LOF_codonPathoRS_CANON",
        assess_forBT_downstreamLOF_CANON : String = TOP_LEVEL_VCF_TAG+"ASSESS_LOF_downstreamPathoRS_CANON",
        
        assess_forBT_geneHasDownstreamLOF : String = TOP_LEVEL_VCF_TAG+"GENEINFO_geneHasPatho_DownstreamLOF",
        assess_forBT_geneHasStartLoss : String = TOP_LEVEL_VCF_TAG+"GENEINFO_geneHasPatho_StartLoss",
        assess_forBT_geneHasStopLoss : String = TOP_LEVEL_VCF_TAG+"GENEINFO_geneHasPatho_StopLoss",
        assess_forBT_geneHasDownstreamLOF_CANON : String = TOP_LEVEL_VCF_TAG+"GENEINFO_geneHasPatho_DownstreamLOF_CANON",
        assess_forBT_geneHasStartLoss_CANON : String = TOP_LEVEL_VCF_TAG+"GENEINFO_geneHasPatho_StartLoss_CANON",
        assess_forBT_geneHasStopLoss_CANON : String = TOP_LEVEL_VCF_TAG+"GENEINFO_geneHasPatho_StopLoss_CANON",
        
        //NEW???????????????
        assess_effectMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_EffectMatch",
        assess_nearEffectMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_NearEffectMatch",
        assess_effectMatchInfo_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_EffectMatch_CANON",
        assess_nearEffectMatchInfo_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_NearEffectMatch_CANON",
        
        assess_pathoEffectMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_EffectMatchRS",
        assess_pathoNearEffectMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_NearEffectMatchRS",
        assess_pathoEffectMatchRS_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_EffectMatchRS_CANON",
        assess_pathoNearEffectMatchRS_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_NearEffectMatchRS_CANON",
        
        assess_benEffectMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvBen_EffectMatchRS",
        assess_benNearEffectMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvBen_NearEffectMatchRS",
        assess_benEffectMatchRS_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvBen_EffectMatchRS_CANON",
        assess_benNearEffectMatchRS_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvBen_NearEffectMatchRS_CANON",
        ////////////
        
        assess_exactMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_ExactMatch",
        assess_aminoMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_AminoMatch",
        assess_nearMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_NearMatch",
        
        assess_pathoExactMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_ExactMatch",
        assess_pathoAminoMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_AminoMatch",
        assess_pathoNearMatchRS : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_NearMatch",
        
        assess_pathoAminoMatchRS_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_AminoMatch_CANON",
        assess_pathoNearMatchRS_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvPath_NearMatch_CANON",
        assess_aminoMatchInfo_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_AminoMatch_CANON",
        assess_nearMatchInfo_CANON : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_NearMatch_CANON",
        
        assess_inSilicoSummary : String = TOP_LEVEL_VCF_TAG+"ACMG_inSilicoSummary",
        
        geneIDs : String = TOP_LEVEL_VCF_TAG+"txGeneIDs",
        
        ec_CallMismatch : String = TOP_LEVEL_VCF_TAG+"EC_CallMismatch",
        ec_CallMismatchStrict : String = TOP_LEVEL_VCF_TAG+"EC_CallMismatch_Strict",
        ec_EnsembleWarnings : String = TOP_LEVEL_VCF_TAG+"EC_ENSEMBLE_WARNINGS",
        ec_alle_callerSets : String = TOP_LEVEL_VCF_TAG+"EC_alle_callerSets",
        ec_singleCallerAllePrefix : String = TOP_LEVEL_VCF_TAG+"EC_altAlle_",
        delims : Vector[String] = Vector(",","|")
      );

  val DefaultVCFAnnoCodes = VCFAnnoCodes();
  
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
    
    def walkVCFFiles(infile : String, outfile : String, chromList : Option[List[String]], infileList : Boolean, vcfCodes : VCFAnnoCodes = VCFAnnoCodes(), verbose : Boolean = true) {
      if(infileList){
        val (peekList,fileiter) = peekIterator(getLinesSmartUnzip(infile),1000);
        val denominator = if(peekList.length == 1000) "?" else peekList.length.toString;
        val (firstIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(peekList.head, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes);
        val vcIter = firstIter ++ fileiter.drop(1).flatMap{file => {
          internalUtils.VcfTool.getVcfIterator(file, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes)._1;
        }}
        val (vcIter2, newHeader) = this.walkVCF(vcIter = vcIter, vcfHeader = vcfHeader)
        val vcfWriter = internalUtils.VcfTool.getVcfWriter(outfile, header = newHeader);
        vcIter2.foreach(vc => {
          vcfWriter.add(vc)
        })
        vcfWriter.close();
      } else {
        walkVCFFile(infile=infile,outfile=outfile,chromList=chromList,vcfCodes=vcfCodes,verbose=verbose);
      }
    }
    
    def chain(walker2 : VCFWalker, flag : Boolean = true) : VCFWalker = {
      if(flag){
        val parent : VCFWalker = this;
        return new VCFWalker {
          def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
            val (iter2,header2) = parent.walkVCF(vcIter = vcIter, vcfHeader = vcfHeader, verbose = verbose);
            walker2.walkVCF(vcIter = iter2, vcfHeader = header2, verbose=verbose);
          }
        }
      } else {
        return this;
      }
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
  
  def splitHeaderLinesByType(headerLines : Seq[VCFHeaderLine]) : (Seq[VCFHeaderLine],Seq[VCFFormatHeaderLine],Seq[VCFInfoHeaderLine]) = {
    val newFormatHeaderLines : (Seq[VCFHeaderLine],Seq[VCFFormatHeaderLine],Seq[VCFInfoHeaderLine]) = headerLines.foldLeft((Seq[VCFHeaderLine](),Seq[VCFFormatHeaderLine](),Seq[VCFInfoHeaderLine]())){ case ((otherLines,fmtLines,infoLines),curr) => {
      curr match {
        case hl : VCFFormatHeaderLine => {
          (otherLines,fmtLines :+ hl, infoLines)
        }
        case hl : VCFInfoHeaderLine => {
          (otherLines,fmtLines,infoLines :+ hl)
        }
        case hl : VCFHeaderLine => {
          (otherLines :+ hl, fmtLines,infoLines)
        }
      }
    }}
    newFormatHeaderLines;
  }
    
  def addHeaderLines(vcfHeader : VCFHeader, headerLines : Seq[VCFHeaderLine]) : VCFHeader = {
    val (newOtherLines,newFmtLines,newInfoLines) = splitHeaderLinesByType(headerLines);
    val (oldOtherLines,oldFmtLines,oldInfoLines) = splitHeaderLinesByType(vcfHeader.getMetaDataInInputOrder().toList);
    
    val oldHeaderLines : List[VCFHeaderLine] = vcfHeader.getMetaDataInInputOrder().toList;
    
    val newFmtSet = newFmtLines.map(hl => hl.getID).toSet;
    val newInfoSet = newInfoLines.map(hl => hl.getID).toSet;
    
    val filtFmtLines = oldFmtLines.filter(hl => ! newFmtSet.contains(hl.getID()));
    val filtInfoLines = oldInfoLines.filter(hl => ! newInfoSet.contains(hl.getID()));
    
    val newHeaderLines : Seq[VCFHeaderLine] = (newOtherLines ++ oldOtherLines ++ filtFmtLines ++ newFmtLines ++ filtInfoLines ++ newInfoLines);
    
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
          if(cs.size == 1){
            val chr = cs.head;
            
            val preProgRep = internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[VariantContext](
                elementTitle  = "skipped lines", lineSec = 60,
                reportFunction  = ((vc : VariantContext, i : Int) => " " + vc.getContig())
            );
            val postProgRep = internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[VariantContext](
                elementTitle  = "processed lines", lineSec = 60,
                reportFunction  = ((vc : VariantContext, i : Int) => " " + vc.getContig())
            );
            
            var i = 1;
            var j = 1;
            vcfIterator.dropWhile{ v => {
              preProgRep.reportProgress(i,v);
              i = i + 1;
              v.getContig() != chr;
            }}.takeWhile{ v => {
              postProgRep.reportProgress(j,v);
              j = j + 1;
              v.getContig() == chr;
            }}
          } else {
            internalUtils.stdUtils.wrapIteratorWithAdvancedProgressReporter[VariantContext](vcfIterator.filter((p : VariantContext) => cs.contains(p.getContig())),
                internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[VariantContext](elementTitle = "lines",lineSec = 60, reportFunction = ((a : VariantContext,i : Int) => a.getContig()))).filter(v =>  cs.contains(v.getContig())); 
          }

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
  
  object SVcfLine {
    
    def readVcfs(lines : Seq[Iterator[String]], withProgress : Boolean = true) : (SVcfHeader,Iterator[SVcfInputVariantLine]) = {
      val (allRawHeaderLines,rawVariantLines0) = lines(0).span(line => line.startsWith("#"));
      val header = readVcfHeader(allRawHeaderLines.toVector);
      
      val rawVariantLines = rawVariantLines0 ++ lines.tail.flatten.filter{line => ! line.startsWith("#")};
      
      val variantLines = if(withProgress){
                         internalUtils.stdUtils.wrapIteratorWithAdvancedProgressReporter[SVcfInputVariantLine](
                           //rawVariantLines.map(line => SVcfInputVariantLine(line,header)),
                           rawVariantLines.map(line => SVcfInputVariantLine(line)), 
                           internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[SVcfInputVariantLine](
                                elementTitle = "lines", lineSec = 60,
                                reportFunction  = ((vc : SVcfInputVariantLine, i : Int) => " " + vc.chrom +" "+ internalUtils.stdUtils.MemoryUtil.memInfo )
                           )
                         )
      } else {
        //rawVariantLines.map(line => SVcfInputVariantLine(line,header))
        rawVariantLines.map(line => SVcfInputVariantLine(line))
      }
      (header,variantLines);
    }
    
    def readVcf(lines : Iterator[String], withProgress : Boolean = true) : (SVcfHeader,Iterator[SVcfInputVariantLine]) = {
      val bufLines = lines.buffered;
      val allRawHeaderLines = extractWhile(bufLines)(line => line.startsWith("#"));
      val header = readVcfHeader(allRawHeaderLines);
      
      val variantLines = if(withProgress){
                         internalUtils.stdUtils.wrapIteratorWithAdvancedProgressReporter[SVcfInputVariantLine](
                           //rawVariantLines.map(line => SVcfInputVariantLine(line,header)),
                           bufLines.map(line => SVcfInputVariantLine(line)), 
                           internalUtils.stdUtils.AdvancedIteratorProgressReporter_ThreeLevelAuto[SVcfInputVariantLine](
                                elementTitle = "lines", lineSec = 60,
                                reportFunction  = ((vc : SVcfInputVariantLine, i : Int) => " " + vc.chrom +" "+ internalUtils.stdUtils.MemoryUtil.memInfo )
                           )
                         )
      } else {
        //rawVariantLines.map(line => SVcfInputVariantLine(line,header))
        bufLines.map(line => SVcfInputVariantLine(line))
      }
      
      (header,variantLines);
    }
    
    //Read Header Lines:
    def readVcfHeader(lines : Seq[String]) : SVcfHeader = {
      var infoLines = Seq[SVcfCompoundHeaderLine]();
      var formatLines = Seq[SVcfCompoundHeaderLine]();
      var otherHeaderLines = Seq[SVcfHeaderLine]();
      
      val (tagLines, nonTagLines) = lines.partition(line => line.startsWith("##"));
      
      
      tagLines.foreach(line => {
        if(line.startsWith("##INFO=")){
          infoLines = infoLines :+ makeCompoundLineFromString(line);
        } else if(line.startsWith("##FORMAT=")){
          formatLines = formatLines :+ makeCompoundLineFromString(line);
        } else {
          otherHeaderLines = otherHeaderLines :+ makeSimpleHeaderLineFromString(line);
        }
      })
      if(nonTagLines.length != 1){
         warning("VCF header line malformed? Found multiple header rows that do not start with \"##\""+ nonTagLines.map{line => "        \""+line+"\""}.mkString("\n"),"MALFORMED_VCF_HEADER_LINE",100);
      }
      
      val titleLine = SVcfTitleLine(nonTagLines.filter{_.startsWith("#CHROM")}.head.split("\t").drop(9));
      
      SVcfHeader(infoLines, formatLines, otherHeaderLines, titleLine);
    }
    
    def makeSimpleHeaderLineFromString(line : String) : SVcfHeaderLine = {
      val tagPair = line.drop(2).split("=",2);
      if(tagPair.length != 2){
        warning("VCF header line malformed?\n\""+line+"\"","MALFORMED_VCF_HEADER_LINE",100);
      }
      val tag = tagPair(0);
      val value = tagPair(1);
      new SVcfHeaderLine(tag,value);
    }
    
    def makeCompoundLineFromString(line : String) : SVcfCompoundHeaderLine = {
      val tagPair = line.drop(2).split("=",2);
      val tag = tagPair(0);
      var hasWarn = false;
      val tagmap = tagPair(1).tail.init.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)").map(compoundString => {
        val ctp = compoundString.split("=",2);
        (ctp(0),(if(ctp.length > 1){
          cleanQuotes(ctp(1))
        }else{
          hasWarn = true;
          warning("Warning: Malformed compound header line:\n\""+line+"\"","Malformed_Compound_Header_Line",-1);
          "."
        }));
      }).toMap;
      
      if(hasWarn){
        warning("SVcfCompoundHeaderLine(tag = "+tag+", ID="+tagmap("ID")+", Number="+tagmap("Number")+", Type="+tagmap("Type")+", desc="+tagmap.getOrElse("Description",".")+")","Malformed_Compound_Header_Line_Continued",-1);
      }
      
      SVcfCompoundHeaderLine(in_tag = tag, ID = tagmap("ID"), Number = tagmap("Number") ,  Type = tagmap("Type") , desc = tagmap.getOrElse("Description","."));
    }
    
    /////////// Read Variant Lines:
    
  }
  
  case class SVcfHeader(var infoLines : Seq[SVcfCompoundHeaderLine], 
                        var formatLines : Seq[SVcfCompoundHeaderLine], 
                        var otherHeaderLines : Seq[SVcfHeaderLine],
                        var titleLine : SVcfTitleLine){
    def getVcfLines : Seq[String] = (otherHeaderLines ++ infoLines ++ formatLines :+ titleLine).map(_.getVcfString);
    
    def addFormatLine(line : SVcfCompoundHeaderLine){
      val idx = formatLines.indexWhere{(p : SVcfCompoundHeaderLine) => {p.ID == line.ID}}
      if( idx != -1 ){
        formatLines = formatLines.updated(idx,line);
      } else {
        formatLines = formatLines :+ line;
      }
    }
    def addInfoLine(line : SVcfCompoundHeaderLine){
      val idx = infoLines.indexWhere{(p : SVcfCompoundHeaderLine) => {p.ID == line.ID}}
      if( idx != -1 ){
        infoLines = infoLines.updated(idx,line);
      } else {
        infoLines = infoLines :+ line;
      }
    }
    
    def getSampleList : Seq[String] = {
      titleLine.sampleList;
    }
    def getSampleIdx(sampID : String) : Int = {
      getSampleList.indexOf(sampID);
    }
    def sampleCt : Int = {
      getSampleList.length;
    }
  }
  
  abstract class SVcfLine {
    def getVcfString : String;
  }
  class SVcfHeaderLine(in_tag : String, var in_value : String) extends SVcfLine{
    var tag : String = in_tag;
    var value : String = in_value;
    
    def getVcfString : String = "##" + tag + "=" + value;
  }
  case class SVcfTitleLine(sampleList : Seq[String]) extends SVcfLine {
    def getVcfString : String = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sampleList.mkString("\t");
  }
  
  case class SVcfCompoundHeaderLine(in_tag : String, var ID : String, var Number : String, var Type : String, var desc : String) extends 
              SVcfHeaderLine(in_tag,"<ID="+ID+",Number="+Number+",Type="+Type+",Description="+addQuotesIfNeeded(desc)+">") {
    //var tag : String = in_tag;
    //var value : String = "<ID="+ID+",Number="+Number+",Type="+Type+",Description=\""+desc+"\">";
  }
  
  /*class SVcfVariantLine extends SVcfLine {
    def getChrom : String;
    def setChrom : Unit;
    def getPos : String;
    def setPos : Unit;
    def getRef : String;
    def setRef : Unit;
    def getAlt : String;
    def setAlt : Unit;
    def 
  }*/
  
  abstract class SVcfVariantLine extends SVcfLine {
    def chrom : String
    def pos : Int
    def id : String
    def ref : String
    def alt : Seq[String]
    def qual : String
    def filter : String
    def info : Map[String,Option[String]];
    def format : Seq[String]
    def genotypes : SVcfGenotypeSet;
    
    def getVcfString : String = chrom + "\t"+pos+"\t"+id+"\t"+ref+"\t"+alt.mkString(",")+"\t"+
                                qual+"\t"+filter+"\t"+info.keySet.toSeq.sorted.map{ case t => {
                                  val v = info(t);
                                  v match {
                                    case Some(sv) => t + "="+sv
                                    case None => t
                                  }
                                }}.mkString(";")+"\t"+
                                format.mkString(":")+"\t"+
                                genotypes.getGenotypeStrings.mkString("\t");
    //def header : SVcfHeader;
    //def getSampleList : Seq[String];

    def getOutputLine() : SVcfOutputVariantLine = {
      SVcfOutputVariantLine(
       in_chrom = chrom,
       in_pos = pos,
       in_id = id,
       in_ref = ref,
       in_alt = alt,
       in_qual = qual,
       in_filter = filter,
       in_info = info,
       in_format = format,
       in_genotypes = genotypes//,
       //in_header = header
      )
    }
    def getSimpleVcfString() : String = chrom + "\t"+pos+"\t"+id+"\t"+ref+"\t"+alt.mkString(",")+"\t"+qual+"\t"+filter;
    
    def getGt(gtTag : String) : Array[String] = {
      error("Op not supported!");
      null;
    }
  }
  
  case class SVcfInputVariantLine(inputLine : String) extends SVcfVariantLine {
    override def getVcfString :  String = inputLine;
    
    lazy val cells = inputLine.split("\t");
    
    lazy val lzy_chrom = cells(0);
    lazy val lzy_pos = string2int(cells(1));
    lazy val lzy_id = cells(2);
    lazy val lzy_ref = cells(3);
    lazy val lzy_alt = cells(4).split(",");
    lazy val lzy_qual = cells(5);
    lazy val lzy_filter = cells(6);
    lazy val lzy_info = cells(7).split(";").map(s => {
      val c = s.split("=",2);
      if(c.length == 1){
        (c(0),None)
      } else {
        (c(0),Some(c(1)));
      }
    }).toMap;
    lazy val lzy_format = cells(8).split(":");
    lazy val lzy_genotypeStrings = cells.drop(9);
    
    def chrom = lzy_chrom;
    def pos = lzy_pos;
    def id = lzy_id;
    def ref = lzy_ref;
    def alt = lzy_alt;
    def qual = lzy_qual;
    def filter = lzy_filter;
    def info = lzy_info;
    def format = lzy_format;
    //def genotypeStrings = lzy_genotypeStrings;
    lazy val lzy_genotypes = SVcfGenotypeSet.getGenotypeSet(lzy_genotypeStrings, format);
    def genotypes = lzy_genotypes;

    override def getGt(gtTag : String) : Array[String] = {
      val gtIdx : Int = format.indexOf(gtTag);
      if(gtIdx == -1) error("GT tag not found!");
      lzy_genotypeStrings.map( gts => (gts.split(":",-1).padTo(gtIdx+1,".")).apply(gtIdx))
    }
    
    //def getSampleList : Seq[String] = in_header.titleLine.sampleList;
    //def header = in_header;
  }
  
  case class SVcfOutputVariantLine(
      var in_chrom : String,
      var in_pos : Int,
      var in_id : String,
      var in_ref : String,
      var in_alt : Seq[String],
      var in_qual : String,
      var in_filter : String,
      var in_info : Map[String,Option[String]],
      in_format : Seq[String],
      var in_genotypes : SVcfGenotypeSet//,
      //var in_header : SVcfHeader
      ) extends SVcfVariantLine {
    def chrom = in_chrom
    def pos = in_pos
    def id = in_id
    def ref = in_ref
    def alt = in_alt
    def qual = in_qual
    def filter = in_filter
    def info = in_info
    def genotypes = in_genotypes
    def format = genotypes.fmt
    //def header = in_header;
    
    def addInfo(infoTag : String, infoVal : String){
      in_info = in_info + ((infoTag,Some(infoVal)));
    }
    
  }
  
  object SVcfGenotypeSet {
    def getGenotypeSet(genotypeStrings : Seq[String], fmt : Seq[String]) : SVcfGenotypeSet = {
      val out = Array.ofDim[String](fmt.length,genotypeStrings.length)
      //val fmtLines = fmt.map(f => {
      //  header.formatLines.find(_.ID == f).get;
      //});
      
      val cells = genotypeStrings.map(_.split(":",-1).padTo(fmt.length,"."));
      
      try{
        cells.indices.foreach(i => {
          (0 until fmt.length).foreach(j => {
            out(j)(i) = cells(i)(j);
          })
        })
      } catch {
        case e : Exception => {
          warning("Error attempting to decode genotypes:\n" + 
                "   FMT: "+fmt.mkString(",") + "\n"+
                "   GENOS: "+genotypeStrings.mkString("\t")+"\n",
                "GENOTYPE_DECODER_ERROR",100
                );
          throw e;
        }
      }
      
      SVcfGenotypeSet(fmt = fmt, genotypeValues = out);
    }
  }
  
  case class SVcfGenotypeSet(var fmt : Seq[String],
                             var genotypeValues : Array[Array[String]]){
    def getGenotypeStrings : Seq[String] = genotypeValues(0).indices.map(j => {
      genotypeValues.indices.map(i => {
        genotypeValues(i)(j)
      }).mkString(":")
    }).toSeq;
    
    def addGenotypeArray(fmtid : String, gval : Array[String]) : Int = {
      val idx = fmt.indexWhere((p : String) => { p == fmtid });
      if(idx != -1){
        genotypeValues(idx) = gval;
        return idx;
      } else {
        fmt = fmt :+ fmtid
        genotypeValues = genotypeValues :+ gval;
        return genotypeValues.length - 1;
      }
    }
    def getPloidy() : Int = genotypeValues(0).map{ (g : String) => {g.split("/").length}}.max;
    
    def addGenotypeArrayIfNew(fmtid : String, gval : Array[String]) : Int = {
      val idx = fmt.indexWhere((p : String) => { p == fmtid });
      if(idx != -1){
        return idx;
      } else {
        fmt = fmt :+ fmtid
        genotypeValues = genotypeValues :+ gval;
        return genotypeValues.length - 1;
      }
    }
    
  }
  

  
  def getSVcfIterator(infile : String, chromList : Option[List[String]],numLinesRead : Option[Int]) : (Iterator[SVcfVariantLine],SVcfHeader) = {
      val (vcfHeader,vcIter) = if(chromList.isEmpty){
        SVcfLine.readVcf(getLinesSmartUnzip(infile),withProgress = true)
      } else if(chromList.get.length == 1){
        val chrom = chromList.get.head;
        SVcfLine.readVcf(getLinesSmartUnzip(infile).filter{line => {
          line.startsWith(chrom+"\t") || line.startsWith("#")
        }},withProgress = true)
      } else {
        val chromSet = chromList.get.toSet;
        val (vh,vi) = SVcfLine.readVcf(getLinesSmartUnzip(infile),withProgress = true)
        (vh,vi.filter(line => { chromSet.contains(line.chrom) }))
      }
      
      val vcIter2 = if(numLinesRead.isDefined){
        vcIter.take(numLinesRead.get);
      } else {
        vcIter
      }
      (vcIter2,vcfHeader);
  }
  def getSVcfIterators(infileString : String, chromList : Option[List[String]],numLinesRead : Option[Int], inputFileList : Boolean = false, withProgress : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader) = {
      val indata = if(inputFileList){
        val (infilePeek,infiles) = peekIterator(getLinesSmartUnzip(infileString),1000);
        val denominator = if(infilePeek.length < 1000) infilePeek.length.toString else "???";
        val headerLines = extractWhile(getLinesSmartUnzip(infilePeek.head).buffered)( a => a.startsWith("#") )
        val allInputLines = flattenIterators(infiles.zipWithIndex.map{case (inf,idx) => addIteratorCloseAction(iter =getLinesSmartUnzip(inf), closeAction = (() => {reportln("finished reading file: "+inf + "("+getDateAndTimeString+")" + "("+(idx+1)+"/"+denominator+")","note")}))}).buffered
        //val headerLines = extractWhile(allInputLines)( a => a.startsWith("#"));
        val remainderLines = allInputLines.filter( a => ! a.startsWith("#"));
        headerLines.iterator ++ remainderLines;
      } else if(infileString.contains(',')){
        val infiles = infileString.split(",");
        val denominator = infiles.length.toString;
        val headerLines = extractWhile(getLinesSmartUnzip(infiles.head).buffered)( a => a.startsWith("#") )
        val allInputLines = flattenIterators(infiles.iterator.zipWithIndex.map{case (inf,idx) => addIteratorCloseAction(iter =getLinesSmartUnzip(inf), closeAction = (() => {reportln("finished reading file: "+inf + "("+getDateAndTimeString+")"+  "("+(idx+1)+"/"+denominator+")","note")}))}).buffered
        //val headerLines = extractWhile(allInputLines)( a => a.startsWith("#"));
        val remainderLines = allInputLines.filter( a => ! a.startsWith("#"));
        headerLines.iterator ++ remainderLines;
      } else {
        getLinesSmartUnzip(infileString)
      }
    
    val (vcfHeader,vcIter) = if(chromList.isEmpty){
        SVcfLine.readVcf(indata,withProgress = withProgress)
      } else if(chromList.get.length == 1){
        val chrom = chromList.get.head;
        SVcfLine.readVcf(indata.filter{line => {
          line.startsWith(chrom+"\t") || line.startsWith("#")
        }},withProgress = withProgress)
      } else {
        val chromSet = chromList.get.toSet;
        val (vh,vi) = SVcfLine.readVcf(indata,withProgress = withProgress)
        (vh,vi.filter(line => { chromSet.contains(line.chrom) }))
      }
      
      val vcIter2 = if(numLinesRead.isDefined){
        vcIter.take(numLinesRead.get);
      } else {
        vcIter
      }
      (vcIter2,vcfHeader);
  }
  
  abstract class SVcfWalker {
    def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader)
    def walkVCFFile(infile :String, outfile : String, chromList : Option[List[String]], numLinesRead : Option[Int]){
      val (vcIter2, vcfHeader) = getSVcfIterator(infile,chromList,numLinesRead);
      
      val (newIter,newHeader) = walkVCF(vcIter2,vcfHeader,verbose=true);
      
      val writer = openWriterSmart(outfile);
      newHeader.getVcfLines.foreach{line => {
        writer.write(line+"\n");
      }}
      newIter.foreach{ line => {
        writer.write(line.getVcfString+"\n");
      }}
      writer.close();
    }
    
    //
    
    def walkVCFFiles(infiles : String, outfile : String, chromList : Option[List[String]], numLinesRead : Option[Int], inputFileList : Boolean){
      val (vcIterRaw, vcfHeader) = getSVcfIterators(infiles,chromList=chromList,numLinesRead=numLinesRead,inputFileList = inputFileList);
      val (newIter,newHeader) = walkVCF(vcIterRaw,vcfHeader,verbose=true);
      val writer = openWriterSmart(outfile);
      newHeader.getVcfLines.foreach{line => {
        writer.write(line+"\n");
      }}
      newIter.foreach{ line => {
        writer.write(line.getVcfString+"\n");
      }}
      writer.close();
    }
  }
  
  
  
  
  
  def parseDebugReport(s : String){
    println(s);
  }
  
  abstract class SFilterLogic[A] {
    def keep(vc : A) : Boolean;
    def printTree() : String;
  }
  
  case class SFilterTrue[A]() extends SFilterLogic[A] {
    def keep(vc : A) : Boolean = true;
    def printTree() : String = "[TRUE]"
  }
  
  case class SFilterByFunc[A](fun : (A => Boolean), params : Seq[String], filterName : String = "UNK") extends SFilterLogic[A] {
    def keep(vc : A) : Boolean = fun(vc);
    def printTree() : String = "["+filterName+" " +params.mkString(" ")+"]";
  }
  
  case class SFilterAND[A](a1: SFilterLogic[A], a2: SFilterLogic[A]) extends SFilterLogic[A] {
    def keep(vc : A) : Boolean = {
      if(! a1.keep(vc)) false;
      else a2.keep(vc);
    }
    def printTree() : String = "["+a1.printTree()+" AND "+a2.printTree()+"]";
  }
  case class SFilterOR[A](a1: SFilterLogic[A], a2: SFilterLogic[A]) extends SFilterLogic[A] {
    def keep(vc : A) : Boolean = {
      if(a1.keep(vc)) true;
      else a2.keep(vc);
    }
    def printTree() : String = "["+a1.printTree()+" OR "+a2.printTree()+"]";
  }
  case class SFilterNOT[A](a1 : SFilterLogic[A]) extends SFilterLogic[A] {
    def keep(vc : A) : Boolean = {
      ! a1.keep(vc);
    }
    def printTree() : String = "[NOT: "+a1.printTree()+"]";
  }

  def findMatchedParenIdx(strs : Seq[String], openParen : String = "(", closeParen : String = ")") : Int = {
      if(strs.head != openParen){
        error("findMatchedParenIdx run on non-paren char!");
      }
      val (depth, idx) = strs.tail.zipWithIndex.iterator.foldLeft((1,1)){ case ((currDepth,currIdx),(currChar,idx)) => {
        if(currDepth == 0){
          (currDepth, currIdx);
        } else if(currChar == openParen){
          (currDepth + 1, currIdx + 1)
        } else if(currChar == closeParen){
          (currDepth - 1, currIdx + 1)
        } else {
          (currDepth, currIdx + 1)
        }
      }}
      idx;
  }
  
  abstract class SFilterLogicParser[A](){
    
    
    
    def parseString(str : String) : SFilterLogic[A] = {
      parseStringArray(str.trim().replaceAll("\\("," \\( ").replaceAll("\\)"," \\) ").trim().split("\\s+").toSeq);
    }
  
    def parseStringArray(strs : Seq[String]) : SFilterLogic[A] = {
      if(strs.length == 0){
        error("Null SVcfFilterLogic!");
      }
      val str = strs.head;
      
      if(str == "("){
        val closeIdx = findMatchedParenIdx(strs);
        val parenStr = strs.take(closeIdx).tail.init;
        parseDebugReport("Entering Parenthetical: [\"" + parenStr.mkString("\",\"") + "\"]" + " with remainder: [\"" + strs.drop(closeIdx).mkString("\",\"")+"\"]");
        parseLogicFunc(parseStringArray(parenStr),strs.drop(closeIdx));
      } else if(str == "NOT"){
        SFilterNOT[A](parseStringArray(strs.tail));
      } else if(strs.head.startsWith("FILT:")) {
        parseDebugReport("Parsing Binary: [" + strs.head +"]");
        parseLogicFunc(parseFilter(strs.head),strs.drop(1));
      } else {
        null;
      }
    }
    
    def parseLogicFunc(a : SFilterLogic[A], strs : Seq[String]) : SFilterLogic[A] = {
      if(strs.length == 0){
        a;
      } else if(strs.head == "AND"){
        SFilterAND[A](a,parseStringArray(strs.tail));
      } else if(strs.head == "OR"){
        SFilterOR[A](a,parseStringArray(strs.tail));
      } else {
        error("Filter logic parse error! Unrecognized logic func: "+strs.head);
        null;
      }
    }
    
    

    def parseFilter(str : String) : SFilterLogic[A] = {
      if(! str.startsWith("FILT:")){
        error("Filter Logic Parse Error! Filters must begin with \"FILT\"");
      }
      val cells = str.split(":")
      if(cells.length < 2){
        error("Filter Logic Parse Error! Filter must have at least 2 elements!");
      }
      val funID = cells(1);
      val params = cells.drop(2)
      val fun = getFilterFunction(cells(1),params);
      SFilterByFunc[A](fun=fun,params=params,filterName=funID);
    }
    
    case class FilterFunction(funcName : String, numParam : Int, desc : String, metaFunc : ((Seq[String]) => ((A) => Boolean))){
      def getFunctionInfo() : String = {
        funcName+"\t"+numParam+"\t"+desc;
      }
    }
    
    def filterFunctionSet : Set[FilterFunction];
    def filterFunctionMap : Map[String,FilterFunction];
    
    def getAllFunctionInfo() : String = {
      var out = "funcName\tnumParam\tdesc\n"
      filterFunctionSet.toSeq.sortBy(f => f.funcName).foreach{ ff => {
        out = out + ff.getFunctionInfo()+"\n";
      }}
      return out;
    }
    
    def getFilterFunction(func : String, params : Seq[String]) : (A => Boolean) = {
      filterFunctionMap.get(func) match {
        case Some(ff) => {
          if(ff.numParam != -1 && params.length != ff.numParam){
            error("Filter Logic Parse Error: Filter Function " + func + " requres "+ff.numParam + " parameters. Found: Params=[\""+params.mkString("\",\"")+"\"]");
          }
          ff.metaFunc(params);
        }
        case None => {
          error("Filter Logic Parse Error: Filter function " + func + " not recognized! Legal functions are: \n"+getAllFunctionInfo());
          return null;
        }
      }
    }
  }
  
  case class SVcfFilterLogicParser() extends SFilterLogicParser[SVcfVariantLine]{
    
    def tagNonMissing(tag : String, a : SVcfVariantLine) : Boolean = {
      a.info.containsKey(tag) && a.info(tag).get != "."
    }
    def tagMissing(tag : String, a : SVcfVariantLine) : Boolean = {
      (! a.info.containsKey(tag)) || a.info(tag).get == "."
    } 
    
    //Map(funcID) = (numParam,desc,metaFunction(params) => function(vc))
    val filterFunctionSetVal : Set[FilterFunction] = Set[FilterFunction]( 
        FilterFunction(funcName="INFO.eq",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && a.info(tag).get == v;
                          }
                        }
                      ), 
        FilterFunction(funcName="INFO.ne",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : SVcfVariantLine) => {
                            (! a.info.containsKey(tag)) || a.info(tag).get != v;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.nm",numParam=1,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a)
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.m",numParam=1,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a)
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.gt",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && string2double(a.info(tag).get) > v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.lt",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && string2double(a.info(tag).get) < v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.ge",numParam=2,desc="", 
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && string2double(a.info(tag).get) >= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.le",numParam=2,desc="", 
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && string2double(a.info(tag).get) <= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.gtm",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a) || string2double(a.info(tag).get) > v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.ltm",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a) || string2double(a.info(tag).get) < v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.gem",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a) || string2double(a.info(tag).get) >= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.lem",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a) || string2double(a.info(tag).get) <= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.in",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && a.info(tag).get.split(",").contains(v) ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.notIn",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a) || a.info(tag).get.split(",").contains(v) ;
                          }
                        }
                      ), 
        FilterFunction(funcName="INFO.inAny", numParam = 2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && a.info(tag).get.split(",").flatMap{s => s.split("\\|")}.contains(v) ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.notInAny", numParam = 2,desc="",
                        (params : Seq[String]) => { 
                          val tag = params(0);
                          val v = params(1);
                          (a : SVcfVariantLine) => {
                            tagMissing(tag,a) || a.info(tag).get.split(",").flatMap{s => s.split("\\|")}.contains(v) ;
                          }
                        }
                      ),
        FilterFunction(funcName="INFO.inAnyOf", numParam = -1,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params.tail.toSet;
                          (a : SVcfVariantLine) => {
                            tagNonMissing(tag,a) && a.info(tag).get.split(",").flatMap{s => s.split("\\|")}.toSet.intersect(v).size > 0 ;
                          }
                        }
                      ),
        FilterFunction(funcName="REF.len.eq", numParam = 1,desc="",
                        (params : Seq[String]) => {
                          val len = string2int(params(0));
                          (a : SVcfVariantLine) => {
                            a.ref.length == len
                          }
                        }
                      ),
        FilterFunction(funcName="ALT.len.eq", numParam = 1,desc="",
                        (params : Seq[String]) => {
                          val len = string2int(params(0));
                          (a : SVcfVariantLine) => {
                            a.alt.head.length == len
                          }
                        }
                      ),
        FilterFunction(funcName="REF.eq", numParam = 1,desc="",
                        (params : Seq[String]) => {
                          val gt = (params(0));
                          (a : SVcfVariantLine) => {
                            a.ref == gt
                          }
                        }
                      ),
        FilterFunction(funcName="ALT.eq", numParam = 1,desc="",
                        (params : Seq[String]) => {
                          val gt = (params(0));
                          (a : SVcfVariantLine) => {
                            a.alt.head == gt
                          }
                        }
                      ),
        FilterFunction(funcName="REF.isOneOf", numParam = -1,desc="",
                        (params : Seq[String]) => {
                          val gtset = params.toSet;
                          (a : SVcfVariantLine) => {
                            gtset.contains(a.ref)
                          }
                        }
                      ),
        FilterFunction(funcName="ALT.isOneOf", numParam = -1,desc="",
                        (params : Seq[String]) => {
                          val gtset = params.toSet;
                          (a : SVcfVariantLine) => {
                            gtset.contains(a.alt.head)
                          }
                        }
                      ),
        FilterFunction(funcName="GENO.hasTagPairMismatch", numParam = 2,desc="",
                        (params : Seq[String]) => {
                          val tag1 = params(0);
                          val tag2 = params(1);
                          (a : SVcfVariantLine) => {
                            val idx1 = a.genotypes.fmt.indexOf(tag1);
                            val idx2 = a.genotypes.fmt.indexOf(tag2);
                            if(idx1 != -1 && idx2 != -1){
                              a.genotypes.genotypeValues(idx1).zip(a.genotypes.genotypeValues(idx2)).exists{ case (v1,v2) => {
                                v1 != v2;
                              }}
                            } else {
                              false;
                            }
                          }
                        }
                      ),
        FilterFunction(funcName="GENO.hasTagPairGtStyleMismatch", numParam = 2,desc="",
                        (params : Seq[String]) => {
                          val tag1 = params(0);
                          val tag2 = params(1);
                          (a : SVcfVariantLine) => {
                            val idx1 = a.genotypes.fmt.indexOf(tag1);
                            val idx2 = a.genotypes.fmt.indexOf(tag2);
                            if(idx1 != -1 && idx2 != -1){
                              a.genotypes.genotypeValues(idx1).zip(a.genotypes.genotypeValues(idx2)).exists{ case (v1,v2) => {
                                (! v1.contains('.')) && (! v2.contains('.')) && v1 != v2;
                              }}
                            } else {
                              false;
                            }
                          }
                        }
                      ),
        FilterFunction(funcName="GENO.MAFgt", numParam = 2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val frac = string2float(params(1));
                          (a : SVcfVariantLine) => {
                            val idx = a.genotypes.fmt.indexOf(tag);
                            if(idx != -1){
                              (a.genotypes.genotypeValues(idx).count{ x => {
                                x.split("/").contains("1")
                              }}.toFloat / a.genotypes.genotypeValues(idx).length.toFloat) > frac
                            } else {
                              false;
                            }
                          }
                        }
                      ),
        FilterFunction(funcName="GENO.MAFlt", numParam = 2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val frac = string2float(params(1));
                          (a : SVcfVariantLine) => {
                            val idx = a.genotypes.fmt.indexOf(tag);
                            if(idx != -1){
                              (a.genotypes.genotypeValues(idx).count{ x => {
                                x.split("/").contains("1")
                              }}.toFloat / a.genotypes.genotypeValues(idx).length.toFloat) < frac
                            } else {
                              false;
                            }
                          }
                        }
                      ),
        FilterFunction(funcName="FILTER.eq",numParam=1,desc="",
                        (params : Seq[String]) => {
                          val v = params(0);
                          (a : SVcfVariantLine) => {
                            a.filter == v;
                          }
                        }
                      ),
        FilterFunction(funcName="FILTER.ne",numParam=1,desc="",
                        (params : Seq[String]) => {
                          val v = params(0);
                          (a : SVcfVariantLine) => {
                            a.filter != v;
                          }
                        }
                      )
    ); 
    
    val filterFunctionMapVal : Map[String,FilterFunction] = {
      filterFunctionSet.map{ ff => {
        (ff.funcName,ff)
      }}.toMap
    }
    def filterFunctionSet : Set[FilterFunction] = filterFunctionSetVal;
    def filterFunctionMap : Map[String,FilterFunction] = filterFunctionMapVal;
  }
  
  case class SGenotypeFilterLogicParser() extends SFilterLogicParser[(SVcfVariantLine,Int)]{
    
    def getTag(tag : String, a : (SVcfVariantLine,Int)) : Option[String] = {
      if(a._1.format.contains(tag)){
        val value = a._1.genotypes.genotypeValues(a._1.format.indexOf(tag))(a._2);
        if(value == ".") None;
        else Some(value);
      } else {
        None;
      }
    }
    def getTagArray(tag : String, a : (SVcfVariantLine,Int)) : Option[String] = {
      if(a._1.format.contains(tag)){
        val value = a._1.genotypes.genotypeValues(a._1.format.indexOf(tag))(a._2);
        if(value == ".") None;
        else Some(value);
      } else {
        None;
      }
    }    
    
    //Map(funcID) = (numParam,desc,metaFunction(params) => function(vc))
    val filterFunctionSetVal : Set[FilterFunction] = Set[FilterFunction]( 
        FilterFunction(funcName="GTAG.eq",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isDefined && value.get == v;
                          }
                        }
                      ), 
        FilterFunction(funcName="GTAG.ne",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = params(1);
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            (! value.isDefined) || value.get != v;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.nm",numParam=1,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          (a : (SVcfVariantLine,Int)) => {
                            getTag(tag,a).isDefined
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.m",numParam=1,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          (a : (SVcfVariantLine,Int)) => {
                            ! getTag(tag,a).isDefined
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.gt",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isDefined && string2double(value.get) > v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.lt",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isDefined && string2double(value.get) < v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.ge",numParam=2,desc="", 
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isDefined && string2double(value.get) >= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.le",numParam=2,desc="", 
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isDefined && string2double(value.get) <= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.gtm",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isEmpty || string2double(value.get) > v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.ltm",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isEmpty || string2double(value.get) < v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.gem",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isEmpty || string2double(value.get) >= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="GTAG.lem",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            value.isEmpty || string2double(value.get) <= v ;
                          }
                        }
                      ),
        FilterFunction(funcName="TAGPAIR.match",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag1 = params(0);
                          val tag2 = params(1);
                          (a : (SVcfVariantLine,Int)) => {
                            val t1 = getTag(tag1,a);
                            val t2 = getTag(tag2,a);
                            t1.isDefined && t2.isDefined && t1 == t2;
                          }
                        }
                      ),
                      
        FilterFunction(funcName="GTAGARRAY.gt",numParam=3,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val idx = string2int(params(1));
                          val v = string2double(params(2));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            if(value.isEmpty) false;
                            else {
                              val arr = value.get.split(",");
                              if(arr.length < idx) false;
                              else if(arr(idx) == ".") false;
                              else {
                                string2double(arr(idx)) > v;
                              }
                            }
                          }
                        }
                      ),
        FilterFunction(funcName="GTAGARRAYSUM.gt",numParam=2,desc="",
                        (params : Seq[String]) => {
                          val tag = params(0);
                          val v = string2double(params(1));
                          (a : (SVcfVariantLine,Int)) => {
                            val value = getTag(tag,a);
                            if(value.isEmpty) false;
                            else {
                              value.get.split(",").map{(s) => if(s == ".") 0.toDouble else string2double(s)}.sum > v
                            }
                          }
                        }
                      )
                      
    ); 
    
    val filterFunctionMapVal : Map[String,FilterFunction] = {
      filterFunctionSet.map{ ff => {
        (ff.funcName,ff)
      }}.toMap
    }
    def filterFunctionSet : Set[FilterFunction] = filterFunctionSetVal;
    def filterFunctionMap : Map[String,FilterFunction] = filterFunctionMapVal;
  }
  
  
  
  /*val vcfOpBinaryFuncMap :  Map[String,((Seq[String]) => ((SVcfVariantLine,String,String) => Boolean))] = Map[String,((Seq[String]) => ((SVcfVariantLine,String,String) => Boolean))](
        ("INFO:eq",(params : Seq[String]) => {
          ((a : SVcfVariantLine, v1 : String, v2 : String) => {
            a.info.containsKey(v1) && a.info(v1).get == v2;
          })
        })
      );*/
  
/*
         ("INFO:ne",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          (! a.info.containsKey(v1)) || a.info(v1).get != v2;
        }),
        ("INFO:nm",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && a.info(v1).get != "."
        }),
        ("INFO:m",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          (! a.info.containsKey(v1)) || a.info(v1).get == "."
        }),
        ("INFO:gt",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && a.info(v1).get != "." && string2double(a.info(v1).get) > string2double(v2);
        }),
        ("INFO:lt",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && a.info(v1).get != "." && string2double(a.info(v1).get) < string2double(v2);
        }),
        ("INFO:gtm",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          (! a.info.containsKey(v1)) || a.info(v1).get == "." || string2double(a.info(v1).get) > string2double(v2);
        }),
        ("INFO:ltm",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && string2double(a.info(v1).get) < string2double(v2);
        }),
        ("INFO:ge",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && a.info(v1).get != "." && string2double(a.info(v1).get) >= string2double(v2);
        }),
        ("INFO:le",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && a.info(v1).get != "." && string2double(a.info(v1).get) <= string2double(v2);
        }),
        ("INFO:gem",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          (! a.info.containsKey(v1)) || a.info(v1).get == "." || string2double(a.info(v1).get) >= string2double(v2);
        }),
        ("INFO:lem",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          (! a.info.containsKey(v1)) || a.info(v1).get == "." || string2double(a.info(v1).get) <= string2double(v2);
        }),
        ("INFO:in",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.info.containsKey(v1) && a.info(v1).get != "." && a.info(v1).get.split(",").contains(v2)
        }),
        ("INFO:notIn",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          (! a.info.containsKey(v1)) || a.info(v1).get == "." || (! a.info(v1).get.split(",").contains(v2))
        }),
        ("FILTER:eq",(a : SVcfVariantLine, v1 : String, v2 : String) => {
          a.filter == v2;
        })
 */
  
  //val SVcfFilterParser = new SFilterLogicParser[SVcfVariantLine](vcfOpBinaryFuncMap);
  
  
}















