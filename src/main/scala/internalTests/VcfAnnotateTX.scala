package internalTests

import htsjdk.variant._;
import htsjdk.variant.variantcontext._;
import htsjdk.variant.vcf._;
import java.io.File;
//import scala.collection.JavaConversions._
import java.io._;
import internalUtils.commandLineUI._;
import internalUtils.Reporter._;

//import scala.collection.JavaConversions._
import scala.collection.JavaConverters._

import internalUtils.optionHolder._;
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import internalUtils.commandLineUI._;
import internalUtils.fileUtils._;
import internalUtils.TXUtil._;
import internalUtils.TXUtil

import internalUtils.VcfTool._;

import htsjdk.variant._;
import htsjdk.variant.variantcontext._;
import htsjdk.variant.vcf._;

import internalUtils.genomicUtils._;
import internalUtils.commonSeqUtils._;

import com.timgroup.iterata.ParIterator.Implicits._;

object VcfAnnotateTX {
  
  class addTXAnno extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "addTxInfoToVCF", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility adds an array of new VCF tags with information about the transcriptional changes caused by each variant. "+BETA_WARNING,
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "summaryFile", 
                                         arg = List("--summaryFile"), 
                                         valueName = "filename.txt",  
                                         argDesc =  "..."
                                        ) ::
                    new UnaryArgument( name = "cdsRegionContainsStop",
                                         arg = List("--cdsRegionContainsStop"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "splitMultiAllelics",
                                         arg = List("--splitMultiAllelics"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::
                                       
                    new BinaryOptionArgument[String](
                                         name = "txInfoFile", 
                                         arg = List("--txInfoFile"), 
                                         valueName = "txInfoFile.txt",  
                                         argDesc =  "Outputs an optional debugging file."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "inputSavedTxFile", 
                                         arg = List("--inputSavedTxFile"), 
                                         valueName = "txdata.data.txt.gz",  
                                         argDesc =  "Loads a saved TXdata file. Either this parameter OR the --genomeFA parameter must be set. Using this file will be much faster than regenerating the tx data from the gtf/fasta."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "outputSavedTxFile", 
                                         arg = List("--outputSavedTxFile"), 
                                         valueName = "txdata.data.txt.gz",  
                                         argDesc =  "Creates a saved TXdata file, for faster loading in future runs. This contains metadata about each transcript in a machine-readable format."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "txToGeneFile", 
                                         arg = List("--txToGeneFile"), 
                                         valueName = "txToGene.txt",  
                                         argDesc =  "File containing the mapping of transcript names to gene symbols. This file must have 2 columns: the txID and the geneID. No header line."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "groupFile", 
                                         arg = List("--groupFile"), 
                                         valueName = "groups.txt",  
                                         argDesc =  "File containing a group decoder. This is a simple 2-column file (no header line). The first column is the sample ID, the 2nd column is the group ID."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "superGroupList", 
                                         arg = List("--superGroupList"), 
                                         valueName = "sup1,grpA,grpB,...;sup2,grpC,grpD,...",  
                                         argDesc =  "A list of top-level supergroups. Requires the --groupFile parameter to be set."
                                        ) :: 
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "genomeFA", 
                                         arg = List("--genomeFA"), 
                                         valueName = "genome.fa.gz",  
                                         argDesc =  "The genome fasta file."
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "input VCF file" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "gtffile.gtf.gz",
                                         argDesc = "A gene annotation GTF file." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output file."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){ 
         runAddTXAnno( parser.get[String]("infile"),
                       parser.get[String]("outfile"),
                       gtffile = parser.get[String]("gtffile"),
                       genomeFA = parser.get[Option[String]]("genomeFA"),
                       txInfoFile = parser.get[Option[String]]("txInfoFile"),
                       summaryFile = parser.get[Option[String]]("summaryFile"),
                       addStopCodon = ! parser.get[Boolean]("cdsRegionContainsStop"),
                       inputSavedTxFile = parser.get[Option[String]]("inputSavedTxFile"),
                       outputSavedTxFile = parser.get[Option[String]]("outputSavedTxFile"),
                       chromList = parser.get[Option[List[String]]]("chromList"),
                       txToGeneFile = parser.get[Option[String]]("txToGeneFile"),
                       groupFile = parser.get[Option[String]]("groupFile"),
                       superGroupList = parser.get[Option[String]]("superGroupList"),
                       splitMultiAllelics = parser.get[Boolean]("splitMultiAllelics")
             )
       }
     }
   } 
  
  
  def runAddTXAnno(vcffile : String, outfile : String, 
                gtffile : String, 
                genomeFA : Option[String],
                summaryFile : Option[String],
                txInfoFile : Option[String],
                addStopCodon : Boolean,
                inputSavedTxFile : Option[String],
                outputSavedTxFile : Option[String],
                chromList : Option[List[String]],
                txToGeneFile : Option[String],
                groupFile : Option[String],
                superGroupList : Option[String],
                splitMultiAllelics : Boolean,
                bufferSize : Int = 32, 
                vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ){
    
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(vcffile, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes);
    
    val summaryWriter = if(summaryFile.isEmpty) None else Some(openWriterSmart(summaryFile.get));
    
    val (vcIter2, newHeader) = 
         AddTxAnnoWalker(
                gtffile =gtffile, 
                genomeFA =genomeFA, 
                summaryWriter =summaryWriter,
                txInfoFile = txInfoFile,
                addStopCodon =addStopCodon,
                inputSavedTxFile =inputSavedTxFile,
                outputSavedTxFile =outputSavedTxFile,
                chromList =chromList,
                txToGeneFile = txToGeneFile,
                bufferSize = bufferSize,
                vcfCodes =vcfCodes
                ).walkVCF(vcIter = vcIter,
                          vcfHeader = vcfHeader)
    
    val (vcIter3,header3) = if(groupFile.isEmpty){
      (vcIter2,newHeader)
    } else {
      AddGroupInfoAnno(groupFile = groupFile, groupList = None, superGroupList  = superGroupList, chromList = chromList).walkVCF(vcIter2,newHeader);
    }
    val (vcIter4,header4) = if(splitMultiAllelics){
      SplitMultiAllelics(vcfCodes = vcfCodes, clinVarVariants = false, splitSimple = false).walkVCF(vcIter3,header3);
    } else {
      (vcIter3,header3);
    }
    
    val vcfWriter = internalUtils.VcfTool.getVcfWriter(outfile, header = header4);
    
    vcIter4.foreach(vc => {
      vcfWriter.add(vc)
    })
    vcfWriter.close();
    
    if(! summaryWriter.isEmpty) summaryWriter.get.close();
  }

  case class FilterChromWalker(chromList : Option[List[String]]) extends internalUtils.VcfTool.VCFWalker {
    val chromSet = chromList match {
      case Some(cl) => {
        Some(cl.toSet);
      }
      case None => {
        None;
      }
    }
    
    val walkerFunc = chromSet match {
      case Some(cs) => {
        (vcIter : Iterator[VariantContext],vcfHeader : VCFHeader) => (vcIter.filter(p => cs.contains(p.getContig())),vcfHeader);
      }
      case None => {
        (vcIter : Iterator[VariantContext],vcfHeader : VCFHeader) => (vcIter,vcfHeader);
      }
    }
    def walkVCF(vcIter : Iterator[VariantContext],vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = walkerFunc(vcIter,vcfHeader);
  }
  
  case class AddTxAnnoWalker(gtffile : String, 
                genomeFA : Option[String], 
                //outfile : String, 
                //summaryFile : Option[String],
                summaryWriter : Option[WriterUtil],
                txInfoFile : Option[String],
                addStopCodon : Boolean,
                inputSavedTxFile : Option[String],
                outputSavedTxFile : Option[String],
                chromList : Option[List[String]],
                txToGeneFile : Option[String],
                bufferSize : Int = 32, 
                vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ) extends internalUtils.VcfTool.VCFWalker {
      
      reportln("Starting TX parse...","progress");
      
      val chromSet = chromList match {
        case Some(lst) => Some(lst.toSet);
        case None => None;
      }
      
      val keepChromFunc = chromSet match {
        case Some(cs) => {
          ((x : String) => cs.contains(x));
        }
        case None => ((x : String) => true);
      }
      
      val TXSeq : scala.collection.mutable.Map[String,TXUtil] = inputSavedTxFile match {
        case Some(txf) => {
          val ipr = internalUtils.stdUtils.IteratorProgressReporter_ThreeLevel("lines", 200, 1000 , 2000 )
          val wrappedIter = internalUtils.stdUtils.wrapIteratorWithProgressReporter(getLinesSmartUnzip(txf) , ipr )
          val txs = scala.collection.mutable.AnyRefMap[String,TXUtil]();
          wrappedIter.foreach{line => {
            val tx = buildTXUtilFromString(line);
            if(keepChromFunc(tx.chrom)) txs += (tx.txID,tx)
          }}
          txs;
        }
        case None => {
          if(genomeFA.isEmpty){
            error("FATAL ERROR: Either the --inputSavedTxFile parameter must be set, OR --genomeFA must be set!")
          }
          buildTXUtilsFromAnnotation(gtffile,genomeFA.get,addStopCodon=addStopCodon,chromSet=chromSet);
        }
      }
      reportln("Finished TX parse.","progress");
      
      outputSavedTxFile match {
        case Some(txf) => {
          val txWriter = openWriterSmart(txf,false);
          TXSeq.foreach{ case (txID,tx) => {
            txWriter.write(tx.saveToString()+"\n");
          }}
          txWriter.close();
        }
        case None => {
          //do nothing!
        }
      }
      
      txInfoFile match {
        case Some(txFile) => {
          val txWriter = openWriterSmart(txFile,false);
          TXSeq.take(100).foreach{ case (txID,tx) => {
            try {
              txWriter.write(tx.toStringVerbose()+"\n");
            } catch {
              case e : Exception => {
                reportln("Caught error on TX:","warn");
                reportln(tx.toStringShort()+"\n","warn");
                throw e;
              }
            }
          }}
          txWriter.close();
        }
        case None => {
          //do nothing
        }
      }
      reportln("Valid Full-Length TX: "+TXSeq.count{case (txID,tx) => {tx.isValidFullLenTX}}+"/"+TXSeq.size,"debug");
      
      reportln("Starting TX GenomicArrayOfSets...","progress");
      val txgaos = new internalUtils.qcGtfAnnotationBuilder(gtffile=gtffile, flatgtffile = None, stranded = false, stdCodes = GtfCodes(), flatCodes = GtfCodes()).txArrayWithBuffer
      reportln("Finished TX GenomicArrayOfSets.","progress");
      
    def walkVCF(vcIter : Iterator[VariantContext],vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      val newHeaderLines = List(
            new VCFInfoHeaderLine(vcfCodes.txList_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "List of known transcripts found to overlap with the variant."),
            new VCFInfoHeaderLine(vcfCodes.vType_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, a |-delimited list indicating the deletion type for each overlapping transcript (see "+vcfCodes.txList_TAG+" for the transcripts and "+vcfCodes.vMutLVL_TAG+" for info on the type description)"),
            new VCFInfoHeaderLine(vcfCodes.vMutG_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, the genomic change, in HGVS format."),
            new VCFInfoHeaderLine(vcfCodes.vMutR_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, a |-delimited list indicating mRNA change (in HGVS format) for each transcript (see "+vcfCodes.txList_TAG+")"),
            new VCFInfoHeaderLine(vcfCodes.vMutC_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, a |-delimited list indicating cDNA change (in HGVS format) for each transcript (see "+vcfCodes.txList_TAG+")"),
            new VCFInfoHeaderLine(vcfCodes.vMutP_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, a |-delimited list indicating amino-acid change (in HGVS format) for each transcript (see "+vcfCodes.txList_TAG+")"),
            new VCFInfoHeaderLine(vcfCodes.vTypeShort_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, the worst amino acid change type found over all transcripts."),
            new VCFInfoHeaderLine(vcfCodes.vMutPShort_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "For each allele, one of the protein changes for the worst variant type found over all transcripts"),
            new VCFInfoHeaderLine(vcfCodes.vMutINFO_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, "Raw variant info for each allele."),
            new VCFInfoHeaderLine(vcfCodes.vMutLVL_TAG, VCFHeaderLineCount.A, VCFHeaderLineType.String, 
                                       "For each allele, the rough deliteriousness level of the variant, over all transcripts. Possible values, in order of deliteriousness "+
                                       "SYNON (synonymous mutation), "+
                                       "PSYNON (Probably-synonymous, indicates that the variant is within a transcript's introns or near a genes endpoints), "+
                                       "UNK (unknown), "+
                                       "NONSYNON (Changes one or more amino acids, or loss of the stop codon), "+
                                       "PLOF (possible loss of function, includes variants that might break splice junctions, and loss of the start codon), and "+
                                       "LLOF (likely loss of function, includes frameshift indels and variants that add early stop codons")
      ) ++ ( if(txToGeneFile.isEmpty) List() else {
        List(
          new VCFInfoHeaderLine(vcfCodes.geneIDs, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene Symbol for each tx.")    
        )
      })
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      
      reportln("Starting VCF read/write...","progress");
      
      val txToGene : Option[(String => String)] = txToGeneFile match {
        case Some(f) => {
          val txToGeneMap = getLinesSmartUnzip(f).map(line => {
            val cells = line.split("\t");
            (cells(0),cells(1));
          }).toMap
          Some(((s : String) => txToGeneMap.getOrElse(s,s)));
        } 
        case None => {
          None;
        }
      }
      
      //val writer = if(summaryFile.isEmpty) None else Some(openWriterSmart(summaryFile.get));
      return ( vcIter.map(v => {
        annotateVcfStream(v,summaryWriter,vcfCodes,bufferSize,txgaos,TXSeq,txToGene)
      }), newHeader );
      
    }
  }
  
  
  def annotateVcfStream(v : VariantContext, writer : Option[WriterUtil], vcfCodes : VCFAnnoCodes,
                        bufferSize : Int, txgaos : GenomicArrayOfSets[String],
                        TXSeq : scala.collection.mutable.Map[String,TXUtil],
                        txToGene : Option[(String => String)]) : VariantContext = {
      var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
      var vTypeList = Vector[Vector[String]]();
      var vMutG = Vector[String]();
      var vMutR = Vector[Vector[String]]();
      var vMutC = Vector[Vector[String]]();
      var vMutP = Vector[Vector[String]]();
      
      var vMutPshort = Vector[String]();
      var vTypeListShort = Vector[String]();
      var vLevel = Vector[String]();
      
      var vInfo = Vector[Vector[String]]();
      
      val refAlle = v.getReference();
      val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a));
      if(altAlleles.length > 0){
        val start = v.getStart - 1
        val end = start + refAlle.length() // math.max(refAlle.length(), altAlleles.map(a => a.length()).max)
        val txList = txgaos.findIntersectingSteps(GenomicInterval(chromName = v.getContig(), strand = '.', start = start - bufferSize, end = end + bufferSize)).foldLeft(Set[String]()){case (soFar,(iv,gset))=>{ soFar ++ gset }}.filter(TXSeq.contains(_)).filter{TXSeq(_).isValidFullLenTX}.toVector.sorted;
        txToGene match {
          case Some(fun) => {
            vb = vb.attribute(vcfCodes.geneIDs, txList.map(x => fun(x)).toList.asJava);
          }
          case None => {
            //do nothing
          }
        }
        if(! writer.isEmpty) {
              writer.get.write(v.getContig()+"\t"+start+"\t"+end+"\t"+""+"\t"+""+"\t"+""+"\t"+""+"\t"+""+"\t"+""+"\t"+""+"\t"+""+"\t \n");
              writer.get.write("\t\t\t"+v.getAttributeAsString("ANN","NA")+"\t"+v.getAttributeAsString("LOF","NA")+"\n");
        }
        
        for(altAlle <- altAlleles){
          val (ref,alt) = (refAlle.getBaseString(), altAlle.getBaseString());
          val mutG = "g."+getMutString(ref, alt, pos = start, getPosString = {(a) => (a + 1).toString},swapStrand=false);
          vMutG = vMutG :+ mutG
          vMutR = vMutR :+ Vector[String]();
          vMutC = vMutC :+ Vector[String]();
          vMutP = vMutP :+ Vector[String]();
          vTypeList = vTypeList :+ Vector[String]();
          vInfo = vInfo :+ Vector[String]();
          
          for(tx <- txList.map(TXSeq(_))){
            val mutR = "r."+tx.getRnaMut(ref,alt,start)  //getMutString(ref,alt, pos = start, getPosString = {(g) => tx.getRPos(g)}, swapStrand = (tx.strand == '-'));
            val mutC = "c."+tx.getCdsMut(ref,alt,start)  //getMutString(ref,alt, pos = start, getPosString = {(g) => tx.getCPos(g)}, swapStrand = (tx.strand == '-'));
            try{
              val mutP = tx.getProteinMut(ref,alt,gPos=start);
              vMutR = vMutR.updated(vMutR.length-1,vMutR.last :+ mutR);
              vMutC = vMutC.updated(vMutC.length-1,vMutC.last :+ mutC);
              
              
              vMutP = vMutP.updated(vMutP.length-1,vMutP.last :+ mutP.pvar);
              vTypeList = vTypeList.updated(vTypeList.length-1,vTypeList.last :+ mutP.varType);
              vInfo = vInfo.updated(vInfo.length-1,vInfo.last :+ mutP.saveToString());
              if(! writer.isEmpty) {
                writer.get.write(v.getContig()+":"+start+"-"+end+"\t"+ref+"\t"+alt+"\t"+tx.txID+"\t"+tx.strand+"\t"+
                                 mutG+"\t"+mutR+"\t"+mutC+"\t"+mutP.pvar+"\t"+mutP.varType+"\t"+mutP.cType+"\t"+mutP.severityType+"\t"+mutP.pType+"\t"+mutP.subType+
                                 "\n");
              }
            } catch {
              case e : Exception => {
                reportln("Error:","warn");
                reportln("TX:","warn");
                reportln(tx.toStringVerbose(),"warn");
                reportln("and variant: ","warn");
                reportln(mutG + "\t"+mutR+"\t"+mutC,"warn");
                throw e;
              }
            }
          }
          
          try{
            val (mutPshort,typeShort,vLvl) = internalUtils.TXUtil.getWorstProteinMut(vMutP.last.zip(vTypeList.last),txList);
            vMutPshort = vMutPshort :+ mutPshort;
            vTypeListShort = vTypeListShort :+ typeShort;
            vLevel = vLevel :+ vLvl;
          } catch {
              case e : Exception => {
                reportln("ERROR:","debug");
                reportln(v.getContig()+":"+start+"-"+end+"\t"+ref+"\t"+alt+"\t"+txList.mkString(",")+"\t"+ vMutP.map(_.mkString("|")).mkString(",")+"\t"+ vTypeList.map(_.mkString("|")).mkString(",")+"\n","debug");
                if( ! writer.isEmpty) writer.get.close();
                throw e;
              }
          }
          
        }
        
        //vb = vb.attribute(vcfCodes.txList_TAG, mkSubDelimString(txList,vcfCodes.delims));
        vb = vb.attribute(vcfCodes.txList_TAG, txList.toList.asJava);
        //vb = vb.attribute(vcfCodes.vType_TAG, mkSubDelimList(vTypeList,vcfCodes.delims).toList.asJava);
        vb = vb.attribute(vcfCodes.vMutG_TAG, vMutG.toList.asJava);
        vb = vb.attribute(vcfCodes.vMutR_TAG, mkSubDelimList(vMutR,vcfCodes.delims));
        vb = vb.attribute(vcfCodes.vMutC_TAG, mkSubDelimList(vMutC,vcfCodes.delims));
        vb = vb.attribute(vcfCodes.vMutP_TAG, mkSubDelimList(vMutP,vcfCodes.delims));
        vb = vb.attribute(vcfCodes.vType_TAG, mkSubDelimList(vTypeList,vcfCodes.delims));
        vb = vb.attribute(vcfCodes.vTypeShort_TAG, vTypeListShort.toList.asJava);
        vb = vb.attribute(vcfCodes.vMutPShort_TAG, vMutPshort.toList.asJava);
        vb = vb.attribute(vcfCodes.vMutLVL_TAG,  vLevel.toList.asJava);
        vb = vb.attribute(vcfCodes.vMutINFO_TAG, mkSubDelimList(vInfo,vcfCodes.delims));
      }
      return vb.make();
  }
  
  class testTXSeqUtil extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "VcfUtilTests", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "Test utility.",   
          argList = 

                    new BinaryArgument[String](name = "subCommand",
                                            arg = List("--cmd"),  
                                            valueName = "cmd", 
                                            argDesc = "The sub-command. This utility serves many separate functions."+
                                                      ""+
                                                      ""+
                                                      "", 
                                            defaultValue = Some("stdtests")
                                       ) :: 
                    new UnaryArgument( name = "cdsRegionContainsStop",
                                         arg = List("--cdsRegionContainsStop"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "version",
                                         arg = List("--version"), // name of value
                                         argDesc = "Flag. If raised, print version." // description
                                       ) ::
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "infile.gtf",
                                         argDesc = "infile" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "genomeFA",
                                         valueName = "genomeFA.fa",
                                         argDesc = "infile" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output file."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
     
     val subCommandList = List[String]("stdtests","extractClinVarLOF");
     
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         testTXSeqUtil(parser.get[String]("gtffile"),
                       parser.get[String]("genomeFA"),
                       parser.get[String]("outfile"),
                       addStopCodon = ! parser.get[Boolean]("cdsRegionContainsStop")
             )
       }
       
     }
   }
  
  def testTXSeqUtil(gtffile : String, genomeFA : String, outfile : String, addStopCodon : Boolean){
    val genomeSeq = internalUtils.genomicAnnoUtils.buildEfficientGenomeSeqContainer(Seq(genomeFA));
    val writer = openWriterSmart(outfile,true);
    
    genomeSeq.shiftBufferTo("chr20",200000);
    
    reportln("genomeSeq test:","debug");
    writer.write("genomeSeq test:\n");
    for(i <- Range(0,100)){
      val (start,end) = (200000 + i * 1000,200000 + (i+1) * 1000);
      writer.write("chr20:"+(start+1)+"-"+end+" "+genomeSeq.getSeqForInterval("chr20",start,end)+"\n");
    }
    //chr20:251503-251908
    writer.write("chr20:"+(251503+1)+"-"+251908+" "+genomeSeq.getSeqForInterval("chr20",251503,251908)+"\n");
    writer.flush();
    //chr20:251504-251908
    
    val TXSeq : scala.collection.mutable.Map[String,TXUtil] = buildTXUtilsFromAnnotation(gtffile,genomeFA, addStopCodon = addStopCodon,debugMode = true);
    
    reportln("Starting output write...","progress");

    
    for(((txID,tx),i) <- TXSeq.iterator.zipWithIndex.take(200)){
      
      val ts = tx.toStringVerbose()
      writer.write(ts+"");

      /*
      writer.write(txID+":\n");
      writer.write(tx.gSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n");
      writer.write(tx.rSpansGS.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n");
      writer.write(tx.rSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n");
      writer.write(tx.seqGS.mkString("")+"\n");
      writer.write(tx.rSpansGS.map{case (i,j) => repString(" ",j-i-1) + "|"}.mkString("")+"\n");
      writer.write(tx.seq.mkString("")+"\n");
      writer.write(tx.rSpans.map{case (i,j) => repString(" ",j-i-1) + "|"}.mkString("")+"\n");
      writer.write(tx.cSeq.mkString("")+"\n");
      */
      /*if(i < 10){
        reportln(".","progress");
        reportln(ts,"progress");
        reportln("ts.length = "+ts.length,"progress");
        //reportln(txID+": "+tx.aa.size,"debug");
        //reportln(tx.aa.mkString("  "),"debug");
      }*/
      //writer.write(internalUtils.commonSeqUtils.getAminoAcidFromSeq(tx.seq,0).mkString("")+"\n");
    }
    reportln("TX anno complete.","progress");
    writer.flush();
    
    val txids = Vector("TESTTX001","TESTTX002","TESTTX003","TESTTX004");
    /*
    val start = 0;
    val end = 140;
    val txList = txids.map(TXSeq(_));
    for(i <- Range(start,end)){
      writer.write(i+"\t"+txList.map(_.getCPos(i)).mkString("\t")+"\n");
      
    }*/
    
    reportln("Output write complete.","progress");
    writer.close();
    reportln("Output writer closed.","progress");
  }
  
  class ConvertGenoPosToCPos extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "ConvertGtoC", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "Test utility.",   
          argList = 
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile.txt",
                                         argDesc = "infile" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "gtffile.gtf",
                                         argDesc = "gtffile" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "genomeFA",
                                         valueName = "genomeFA.fa",
                                         argDesc = "genomeFA" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output file."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
     
     val subCommandList = List[String]("stdtests","extractClinVarLOF");
     
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         ConvertGenoPosToCPos(parser.get[String]("infile"),
                       parser.get[String]("gtffile"),
                       parser.get[String]("genomeFA"),
                       parser.get[String]("outfile")
             )
       }
     }
   }
  
  def ConvertGenoPosToCPos(infile : String, gtffile : String, genomeFA : String,  outfile : String){
    val TXSeq : scala.collection.mutable.Map[String,TXUtil] = buildTXUtilsFromAnnotation(gtffile,genomeFA);
    val writer = openWriterSmart(outfile,true);
      
    //
    
    val lines = getLinesSmartUnzip(infile,true);
    
    for(line <- lines){
      val cells = line.split("\\s+");
      val txid = cells(0);
      val pos = string2int(cells(1));
      writer.write(line+"\t"+TXSeq(txid).getCPos(pos));
    }
    
    writer.close();
  }
  
  

  class ConvertAminoRangeToGenoRange extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "ConvertAminoRangeToGenoRange", 
          quickSynopsis = "", 
          synopsis = "", 
          description = " " + ALPHA_WARNING,   
          argList = 
                    new BinaryArgument[Int](   name = "txCol",
                                               arg = List("--txCol"),  
                                               valueName = "0", 
                                               argDesc = "", 
                                               defaultValue = Some(0)
                                           ) ::
                    new BinaryArgument[Int](   name = "startCol",
                                               arg = List("--startCol"),  
                                               valueName = "1", 
                                               argDesc = "", 
                                               defaultValue = Some(1)
                                           ) ::
                    new BinaryArgument[Int](   name = "endCol",
                                               arg = List("--endCol"),  
                                               valueName = "2", 
                                               argDesc = "", 
                                               defaultValue = Some(2)
                                           ) ::
                    new BinaryOptionArgument[String](
                                         name = "badSpanFile", 
                                         arg = List("--badSpanFile"), 
                                         valueName = "file.txt",  
                                         argDesc =  "..."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "unkTxFile", 
                                         arg = List("--unkTxFile"), 
                                         valueName = "file.txt",  
                                         argDesc =  "..."
                                        ) ::     
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile.txt",
                                         argDesc = "infile" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "txdataFile",
                                         valueName = "txdataFile.txt.gz",
                                         argDesc = "txdataFile" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output file."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
     
     val subCommandList = List[String]("stdtests","extractClinVarLOF");
     
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         ConvertAminoRangeToGenoRange(parser.get[String]("infile"),
                       parser.get[String]("txdataFile"),
                       parser.get[String]("outfile"),
                       parser.get[Int]("txCol"),
                       parser.get[Int]("startCol"),
                       parser.get[Int]("endCol"),
                       badSpanFile =  parser.get[Option[String]]("badSpanFile"),
                       unkTxFile =  parser.get[Option[String]]("unkTxFile")
             )
       }
     }
   }
  
  
  def ConvertAminoRangeToGenoRange(infile : String, txdataFile : String, outfile : String,
                                   txIdx : Int = 0, startIdx : Int = 1, endIdx : Int = 2, 
                                   badSpanFile : Option[String] = None,
                                   unkTxFile : Option[String] = None){
    //val TXSeq : scala.collection.mutable.Map[String,TXUtil] = buildTXUtilsFromAnnotation(gtffile,genomeFA);
    
    val TXSeq : scala.collection.mutable.Map[String,TXUtil] = scala.collection.mutable.AnyRefMap[String,TXUtil]();
    
    reportln("Reading TX file...","progress");
    wrapIteratorWithAdvancedProgressReporter(getLinesSmartUnzip(txdataFile),AdvancedIteratorProgressReporter_ThreeLevelAuto[String]()).foreach(line => {
      val tx = buildTXUtilFromString(line);
      TXSeq.put(tx.txID,tx);
    });
    reportln("Finished with TX file.","progress");
      
    val writer = openWriterSmart(outfile,true);

    val (firstLine,lines) = peekIterator(wrapIteratorWithAdvancedProgressReporter(getLinesSmartUnzip(infile,true),AdvancedIteratorProgressReporter_ThreeLevelAuto[String]()));
    val firstCells = firstLine.split("\t");
    val copyCols = (Range(0,firstCells.size).toSet -- Set(txIdx,startIdx,endIdx)).toVector.sorted;
    
    writer.write("#chrom\tstart\tend\tprotSymbol\tprotLen\ttxID\tdomainStartAA\tdomainEndAA\tdomainID\ttxLenAA\ttxStart\ttxEnd\n")
    
    var badSpanSet = Set[(String,Int,Int,String)]();
    var unkTx = Set[String]();
    for(line <- lines){
      val cells = line.split("\t");
      val txid  = cells(txIdx);
      val protLenAA = string2int(cells(1));
      //val txLenAA = string2int(cells(7));
      val startAA = string2int(cells(startIdx)) - 1;
      val endAA   = if(cells(endIdx).head == '>') string2int(cells(endIdx).tail) else string2int(cells(endIdx))
      val startC = startAA * 3;
      val endC   = endAA * 3;
      
      if(TXSeq.contains(txid)){
          TXSeq.get(txid) match {
            case Some(tx) => {
              val txLenAA = tx.cLen/3;
              if(protLenAA + 1 == txLenAA){
                val gSpans = tx.convertCSpanToGSpan(startC,endC);
              
                gSpans.foreach{case (s,e) => {
                  writer.write(tx.chrom+"\t"+s+"\t"+e+"\t"+line+"\t"+txLenAA+"\t"+tx.gStart+"\t"+tx.gEnd+"\n");
                }}
              } else {
                badSpanSet = badSpanSet + ((tx.chrom,tx.gStart,tx.gEnd,txid)); 
              }
            }
            case None => {
              writer.write(line+"\t-1\t-1\t-1\n");
            }
          }
      } else {
        unkTx = unkTx + txid;
      }
      //writer.write(line+"\t"+TXSeq(txid).getCPos(pos));
    }
    badSpanFile match {
      case Some(f) => {
        val w = openWriterSmart(f);
        w.write("chrom\tstart\tend\ttxID")
        val badSpanList = badSpanSet.toVector.sorted;
        badSpanList.foreach{case (chrom,start,end,txid) => {
          w.write(chrom+"\t"+start+"\t"+end+"\t"+txid+"\n");
        }}
        w.close();
      }
      case None => {
        //do nothing
      }
    }
    unkTxFile match {
      case Some(f) => {
        val w = openWriterSmart(f);
        //w.write("chrom\tstart\tend\ttxID")
        val txList = unkTx.toVector.sorted;
        txList.foreach{txid => {
          w.write(txid+"\n");
        }}
        w.close();
      }
      case None => {
        //do nothing
      }
    }
    
    writer.close();
  }
  
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  class AddGroupSummaries extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "addGroupSummaries", 
          quickSynopsis = "", 
          synopsis = "", 
          description = " " + ALPHA_WARNING,
          argList = 
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "..."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "groupFile", 
                                         arg = List("--groupFile"), 
                                         valueName = "file.txt",  
                                         argDesc =  "..."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "groupList", 
                                         arg = List("--groupList"), 
                                         valueName = "grpA,A1,A2,...;grpB,B1,...",  
                                         argDesc =  "..."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "superGroupList", 
                                         arg = List("--superGroupList"), 
                                         valueName = "sup1,grpA,grpB,...;sup2,grpC,grpD,...",  
                                         argDesc =  "..."
                                        ) ::
                    new UnaryArgument( name = "noCts",
                                         arg = List("--noCts"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) :: 
                    new UnaryArgument( name = "noFrq",
                                         arg = List("--noFrq"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) :: 
                    new UnaryArgument( name = "noAlle",
                                         arg = List("--noAlle"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::     
                    new UnaryArgument( name = "noGeno",
                                         arg = List("--noGeno"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) :: 
                    new UnaryArgument( name = "noMiss",
                                         arg = List("--noMiss"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) :: 
                    new FinalArgument[String](
                                         name = "invcf",
                                         valueName = "variants.vcf",
                                         argDesc = "infput VCF file" // description
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
       AddGroupInfoAnno(groupFile = parser.get[Option[String]]("groupFile"),
                            groupList = parser.get[Option[String]]("groupList"),
                            superGroupList = parser.get[Option[String]]("superGroupList"),
                            chromList = parser.get[Option[List[String]]]("chromList"),
                            addCounts = ! parser.get[Boolean]("noCts"),
                            addFreq   = ! parser.get[Boolean]("noFrq"),
                            addMiss   = ! parser.get[Boolean]("noMiss"),
                            addAlle   = ! parser.get[Boolean]("noAlle"),
                            addHetHom = ! parser.get[Boolean]("noGeno")).walkVCFFile(
                            parser.get[String]("invcf"), 
                            parser.get[String]("outvcf"),
                            chromList = parser.get[Option[List[String]]]("chromList")
                            );
       
       /*runAddGroupInfoAnno( parser.get[String]("invcf"),
                            parser.get[String]("outvcf"),
                            groupFile = parser.get[Option[String]]("groupFile"),
                            groupList = parser.get[Option[String]]("groupList"),
                            superGroupList = parser.get[Option[String]]("superGroupList"),
                            chromList = parser.get[Option[List[String]]]("chromList"),
                            addCounts = ! parser.get[Boolean]("noCts"),
                            addFreq   = ! parser.get[Boolean]("noFrq"),
                            addMiss   = ! parser.get[Boolean]("noMiss"),
                            addAlle   = ! parser.get[Boolean]("noAlle"),
                            addHetHom = ! parser.get[Boolean]("noGeno")                            
           )*/
     }
  }
  }
  
  /*def runAddGroupInfoAnno(infile : String, outfile : String, groupFile : Option[String], groupList : Option[String], superGroupList  : Option[String],
                          chromList : Option[List[String]], 
                          addCounts : Boolean = true, addFreq : Boolean = true, addMiss : Boolean = true, 
                          addAlle : Boolean= true, addHetHom : Boolean = true, 
                          sepRef : Boolean = true, countMissing : Boolean = true,
                          vcfCodes : VCFAnnoCodes = VCFAnnoCodes()){
     
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(infile, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes);
        
    val (vcIter2, newHeader) = AddGroupInfoAnno(groupFile=groupFile,groupList=groupList,superGroupList=superGroupList,chromList=chromList,
                                                addCounts = addCounts, addFreq=addFreq, addMiss=addMiss,
                                                addAlle=addAlle, addHetHom=addHetHom, sepRef=sepRef,
                                                vcfCodes=vcfCodes).walkVCF(vcIter,vcfHeader);
    
    val vcfWriter = internalUtils.VcfTool.getVcfWriter(outfile, header = newHeader);
    
    vcIter2.foreach(vc => {
      vcfWriter.add(vc)
    })
    vcfWriter.close();
  }*/
  
  case class AddGroupInfoAnno(groupFile : Option[String], groupList : Option[String], superGroupList  : Option[String], chromList : Option[List[String]], 
             addCounts : Boolean = true, addFreq : Boolean = true, addMiss : Boolean = true, 
             addAlle : Boolean= true, addHetHom : Boolean = true, 
             sepRef : Boolean = true, countMissing : Boolean = true,
             vcfCodes : VCFAnnoCodes = VCFAnnoCodes()) extends internalUtils.VcfTool.VCFWalker {
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext], VCFHeader) = {
  
      val sampleToGroupMap = new scala.collection.mutable.AnyRefMap[String,Set[String]](x => Set[String]());
      val groupToSampleMap = new scala.collection.mutable.AnyRefMap[String,Set[String]](x => Set[String]());
      var groupSet : Set[String] = Set[String]();
      
      groupFile match {
        case Some(gf) => {
          val r = getLinesSmartUnzip(gf).drop(1);
          r.foreach(line => {
            val cells = line.split("\t");
            if(cells.length != 2) error("ERROR: group file must have exactly 2 columns. sample.ID and group.ID!");
            groupToSampleMap(cells(1)) = groupToSampleMap(cells(1)) + cells(0);
            sampleToGroupMap(cells(0)) = sampleToGroupMap(cells(0)) + cells(1);
            groupSet = groupSet + cells(1);
          })
        }
        case None => {
          //do nothing
        }
      }
      
      groupList match {
        case Some(g) => {
          val r = g.split(";");
          r.foreach(grp => {
            val cells = grp.split(",");
            val grpID = cells.head;
            cells.tail.foreach(samp => {
              sampleToGroupMap(samp) = sampleToGroupMap(samp) + grpID;
              groupToSampleMap(grpID) = groupToSampleMap(grpID) + samp;
            })
            groupSet = groupSet + grpID;
          })
        }
        case None => {
          //do nothing
        }
      }
      superGroupList match {
        case Some(g) => {
          val r = g.split(";");
          r.foreach(grp => {
            val cells = grp.split(",");
            val grpID = cells.head;
            cells.tail.foreach(subGrp => {
              groupToSampleMap(grpID) = groupToSampleMap(grpID) ++ groupToSampleMap(subGrp);
            })
            groupSet = groupSet + grpID;
          })
        }
        case None => {
          //do nothing
        }
      }
  
      val sampNames = vcfHeader.getSampleNamesInOrder().asScala.toVector;
      val groups = groupSet.toVector.sorted;
  
      reportln("Final Groups:","debug");
      for((g,i) <- groups.zipWithIndex){
        reportln("Group "+g+" ("+groupToSampleMap(g).size+")","debug");
      }
      
      //addCounts : Boolean = true, addFreq : Boolean = true, addAlle, addHetHom : Boolean = true, sepRef : Boolean = true,
      val forEachString = "for each possible allele (including the ref)"; //if(sepRef) "for each possible ALT allele (NOT including the ref)" else "for each possible allele (including the ref)";
      val countingString = if(countMissing) " counting uncalled alleles." else " not counting uncalled alleles."
      val newHeaderLines = groups.map(g => {
        (if(addCounts && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpAC_TAG + g, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of alleles called "+forEachString+","+countingString)) else List[VCFHeaderLine]()) ++
        (if(addFreq   && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpAF_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the proportion of alleles called "+"for each alt allele (not including ref)"+","+countingString)) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom) List(new VCFInfoHeaderLine(vcfCodes.grpHomCt_TAG + g, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of homozygous genotypes called "+forEachString+","+countingString   )) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom) List(new VCFInfoHeaderLine(vcfCodes.grpHetCt_TAG + g, VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of heterozygous genotypes called "+forEachString+","+countingString )) else List[VCFHeaderLine]()) ++
        (if(addFreq   && addHetHom) List(new VCFInfoHeaderLine(vcfCodes.grpHomFrq_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the fraction of homozygous genotypes called "+"for each alt allele (not including ref)"+","+countingString  )) else List[VCFHeaderLine]()) ++
        (if(addFreq   && addHetHom) List(new VCFInfoHeaderLine(vcfCodes.grpHetFrq_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the fraction of heterozygous genotypes called "+"for each alt allele (not including ref)"+","+countingString)) else List[VCFHeaderLine]()) ++
        (if(addCounts && addMiss)   List(new VCFInfoHeaderLine(vcfCodes.grpMisCt_TAG + g,  1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of alleles missing (non-called).")) else List[VCFHeaderLine]()) ++
        (if(addFreq   && addMiss)   List(new VCFInfoHeaderLine(vcfCodes.grpMisFrq_TAG + g, 1, VCFHeaderLineType.Float, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the fraction of alleles missing (non-called).")) else List[VCFHeaderLine]()) ++
        List[VCFHeaderLine]()
      }).flatten.toList
  
        //      grpMisCt_TAG : String = TOP_LEVEL_VCF_TAG+"MisCt_GRP_",
        //  grpMisFrq_TAG : String = TOP_LEVEL_VCF_TAG+"MisFrq_GRP_",
        //  grpRefAC_TAG : String = TOP_LEVEL_VCF_TAG+"RefAC_GRP_",
        //  grpRefAF_TAG : String = TOP_LEVEL_VCF_TAG+"RefAF_GRP_",
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      
      return (vcIter.map(vc => {
        var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(vc);
        val gt = vc.getGenotypes().iterateInSampleNameOrder(sampNames.asJava).asScala.toVector.zip(sampNames);
        val alleles = internalUtils.VcfTool.getAllelesInOrder(vc).toVector;
        
        val groupAlleCts = alleles.map(alle => { groups.map(grp => {
            val sampSet = groupToSampleMap(grp);
            gt.foldLeft(0){case (soFar,(g,samp)) => {
              if(sampSet.contains(samp)){
                soFar + g.countAllele(alle);
              } else {
                soFar;
              }
            }}
          })
        });
        
        val groupMisCts = groups.map(grp => {
            val sampSet = groupToSampleMap(grp);
            gt.foldLeft(0){case (soFar,(g,samp)) => {
              if(sampSet.contains(samp)){
                soFar + g.countAllele(Allele.NO_CALL);
              } else {
                soFar;
              }
            }}
          })
        
        val groupHomCts = alleles.map(alle => { groups.map(grp => {
            val sampSet = groupToSampleMap(grp);
            gt.foldLeft(0){case (soFar,(g,samp)) => {
              if(sampSet.contains(samp) && g.countAllele(alle) == 2){
                soFar + 1;
              } else {
                soFar;
              }
            }}
          })
        });
        val groupHetCts = alleles.map(alle => { groups.map(grp => {
            val sampSet = groupToSampleMap(grp);
            gt.foldLeft(0){case (soFar,(g,samp)) => {
              if(sampSet.contains(samp) && g.countAllele(alle) > 1){
                soFar + 1;
              } else {
                soFar;
              }
            }}
          })
        });
        
  
        val groupGenoSums = groups.map(grp => groupToSampleMap(grp).size.toDouble )
        val groupAlleAllSums = groups.zipWithIndex.map{case (grp,gi) => {
            val sampSet = groupToSampleMap(grp);
            gt.foldLeft(0){case (soFar,(g,samp)) => {
              if(sampSet.contains(samp)){
                soFar + g.getPloidy();
              } else {
                soFar;
              }
            }}.toDouble
        }}
        val groupAlleSums = if(countMissing){
          groupAlleAllSums;
        } else {
          groups.zipWithIndex.map{case (grp,gi) => {
            groupAlleCts.map(_(gi)).sum.toDouble
          }}
        }
        
        val groupAlleAF = groupAlleCts.map{groupAC => {
          groupAC.zip(groupAlleSums).map{case (ct,sumCt) => ct.toDouble / sumCt};
        }}
        
        val groupHomFrq = groupHomCts.map{groupCt => {
          groupCt.zip(groupGenoSums).map{case (ct,sumCt) => { ct.toDouble / sumCt }}
          //groupCt.map(_.toDouble / sumCt);
        }}
        val groupHetFrq = groupHetCts.map{groupCt => {
          groupCt.zip(groupGenoSums).map{case (ct,sumCt) => { ct.toDouble / sumCt }}
        }}
        val groupMisFrq = groupMisCts.zip(groupAlleAllSums).map{case (ct,sumCt) => { ct.toDouble / sumCt }}
        
        for((g,i) <- groups.zipWithIndex){
          if(addCounts && addAlle) vb = vb.attribute(vcfCodes.grpAC_TAG + g, groupAlleCts.map(_(i)).toList.asJava);
          if(addFreq && addAlle) vb = vb.attribute(vcfCodes.grpAF_TAG + g, groupAlleAF.tail.map(_(i)).toList.asJava);
          if(addCounts && addHetHom) vb = vb.attribute(vcfCodes.grpHomCt_TAG  + g, groupHomCts.map(_(i)).toList.asJava);
          if(addCounts && addHetHom) vb = vb.attribute(vcfCodes.grpHetCt_TAG  + g, groupHetCts.map(_(i)).toList.asJava);
          if(addFreq && addHetHom) vb = vb.attribute(vcfCodes.grpHomFrq_TAG + g, groupHomFrq.tail.map(_(i)).toList.asJava);
          if(addFreq && addHetHom) vb = vb.attribute(vcfCodes.grpHetFrq_TAG + g, groupHetFrq.tail.map(_(i)).toList.asJava);
          if(addCounts   && addMiss) vb = vb.attribute(vcfCodes.grpMisCt_TAG + g,  groupMisCts(i));
          if(addFreq     && addMiss) vb = vb.attribute(vcfCodes.grpMisFrq_TAG + g, groupMisFrq(i));
        }
        
        vb.make();
      }),newHeader);
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

  class CmdSplitMultiAllelics extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "splitMultiAllelics", 
          quickSynopsis = "", 
          synopsis = "", 
          description = " " + ALPHA_WARNING,   
          argList = 
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "..."
                                        ) ::
                    new UnaryArgument( name = "clinVarVariants",
                                         arg = List("--clinVarVariants"), // name of value
                                         argDesc = ""+
                                                   ""
                                       ) ::
                    new FinalArgument[String](
                                         name = "invcf",
                                         valueName = "variants.vcf",
                                         argDesc = "infput VCF file" // description
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
       runSplitMultiAllelics( parser.get[String]("invcf"),
                              parser.get[String]("outvcf"),
                              chromList = parser.get[Option[List[String]]]("chromList"),
                              clinVarVariants = parser.get[Boolean]("clinVarVariants")
           )
     }
  }
  }
  
  def runSplitMultiAllelics(infile : String, outfile : String, chromList : Option[List[String]], clinVarVariants : Boolean, vcfCodes : VCFAnnoCodes = VCFAnnoCodes()){
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(infile, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes);
        
    val (vcIter2, newHeader) = SplitMultiAllelics(vcfCodes = vcfCodes, clinVarVariants = clinVarVariants, splitSimple = clinVarVariants).walkVCF(vcIter,vcfHeader);
    
    val vcfWriter = internalUtils.VcfTool.getVcfWriter(outfile, header = newHeader);
    
    vcIter2.foreach(vc => {
      vcfWriter.add(vc)
    })
    vcfWriter.close();
  }
  
  val multiAllelicIndexStream = Iterator.from(0).toStream.map(i => {
    Range(0,i).flatMap(k => Range(k,i).map(z => (k,z)))
  }) 
  
  case class SplitMultiAllelics(vcfCodes : VCFAnnoCodes = VCFAnnoCodes(), clinVarVariants : Boolean, splitSimple : Boolean = false) extends internalUtils.VcfTool.VCFWalker {
      val infoCLN = if(clinVarVariants){
        Set[String]("CLNHGVS","CLNSRC","CLNORIGIN","CLNSRCID","CLNSIG","CLNDSDB","CLNDSDBID","CLNREVSTAT","CLNACC","CLNDBN");
      } else {
        Set[String]();
      }
      
    val alleCtCt = new scala.collection.mutable.HashMap[Int, Int]().withDefaultValue(0)
    val badCtCt = new scala.collection.mutable.HashMap[Int, Int]().withDefaultValue(0)
      
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext], VCFHeader) = {
      val newHeaderLines = List[VCFHeaderLine](
          //new VCFInfoHeaderLine(vcfCodes.isSplitMulti_TAG, 0, VCFHeaderLineType.Flag,    "Indicates that this line was split apart from a multiallelic VCF line."),
          new VCFInfoHeaderLine(vcfCodes.splitIdx_TAG,     1, VCFHeaderLineType.Integer, "Indicates the index of this biallelic line in the set of biallelic lines extracted from the same multiallelic VCF line."),
          new VCFInfoHeaderLine(vcfCodes.numSplit_TAG,     1, VCFHeaderLineType.Integer, "Indicates the number of biallelic lines extracted from the same multiallelic VCF line as this line."),
          new VCFInfoHeaderLine(vcfCodes.splitAlle_TAG,     VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "The original set of alternative alleles.")
      ) ++ (if(clinVarVariants){
        List[VCFHeaderLine](
            new VCFInfoHeaderLine("CLNPROBLEM",     1, VCFHeaderLineType.Integer, "Indicates whether the splitting of the Clinvar CLN tags was successful. 1=yes, 0=no.")
        ) ++ infoCLN.map(key => {
            new VCFInfoHeaderLine("ORIG"+key,     VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "")
        })
      } else {List[VCFHeaderLine]()})
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      
      val infoA = newHeader.getMetaDataInInputOrder().asScala.filter(_ match {
        case x : VCFInfoHeaderLine => {
          x.getCountType() == VCFHeaderLineCount.A;
        }
        case _ => {
          false;
        }
      }).map(_ match {
        case x : VCFInfoHeaderLine => {
          x.getID()
        }
      })
      val infoR = newHeader.getMetaDataInInputOrder().asScala.filter(_ match {
        case x : VCFInfoHeaderLine => {
          x.getCountType() == VCFHeaderLineCount.R;
        }
        case _ => {
          false;
        }
      }).map(_ match {
        case x : VCFInfoHeaderLine => {
          x.getID()
        }
      })


      infoCLN.foreach(x => {
        reportln("INFO Line "+x+" is of type CLN","debug");
      })
      infoA.foreach(x => {
        reportln("INFO Line "+x+" is of type A","debug");
      })
      infoR.foreach(x => {
        reportln("INFO Line "+x+" is of type R","debug");
      })
      
      val sampNames = vcfHeader.getSampleNamesInOrder();
      return (addIteratorCloseAction(vcIter.flatMap(vc => {
        if(vc.getNAlleles() <= 2 && (! splitSimple)){
          var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(vc);
          vb = vb.attribute(vcfCodes.numSplit_TAG, "1" );
          vb = vb.attribute(vcfCodes.splitIdx_TAG, "0" );
          alleCtCt(1) = alleCtCt(1) + 1;
          List[VariantContext](vb.make());
        } else {
          val alleles = internalUtils.VcfTool.getAllelesInOrder(vc).toVector;
          val numAlle = alleles.length;
          val refAlle = alleles.head;
          val altAlleles = alleles.tail;
          alleCtCt(altAlleles.length) = alleCtCt(altAlleles.length) + 1;
          val attr = vc.getAttributes().asScala.map{case (key,obj) => { (key,vc.getAttributeAsList(key).asScala.map(_.toString).mkString(",")  ) }}.toMap
          val (copyAttrA,tempAttrSet)    = attr.partition{case (key,obj) => { infoA.contains(key) }};
          val (copyAttrCLN,tempAttrSet2)    = tempAttrSet.partition{case (key,obj) => { infoCLN.contains(key) }};
          val (copyAttrR,simpleCopyAttr) = tempAttrSet2.partition{case (key,obj) => { infoR.contains(key) }};
          
          val attrA = copyAttrA.map{case (key,attr) => {
            val a = {
              val atr = attr.split(",");
              if(atr.length != numAlle - 1){
                warning("WARNING: ATTR: \""+key+"\"=\""+attr+"\" of type \"A\"\n   VAR:\n      atr.length() = "+atr.length + " numAlle = "+numAlle+"\n   "+vc.toStringWithoutGenotypes(),"INFO_LENGTH_WRONG",10);
                repToSeq(atr.mkString(","),numAlle - 1);
              } else {
                atr.toSeq;
              }
            }
            (key,a)
          }}
          val attrR = copyAttrR.map{case (key,attr) => {
            val a = {
              val atr = attr.split(",");
              if(atr.length != numAlle){
                warning("WARNING: ATTR: \""+key+"\"=\""+attr+"\" of type \"R\"\n   VAR:\n      atr.length() = "+atr.length + " numAlle = "+numAlle+"\n   "+vc.toStringWithoutGenotypes(),"INFO_LENGTH_WRONG",10);
                repToSeq(atr.mkString(","),numAlle);
              } else {
                atr.toSeq;
              }
            }
            (key,a)
          }}

          val ANNSTR =  trimBrackets(attr.getOrElse("ANN","")).split(",").map(ann => { (ann.trim.split("\\|")(0),ann.trim) });
          
          //warning("ANNSTR: \n      "+ANNSTR.map{case (a,b) => { "(\""+a+"\",\""+b+"\")" }}.mkString("\n      ")+"","ANN_FIELD_TEST",100);
          
          val ANN = if(attr.contains("ANN")) Some(altAlleles.map(_.getBaseString()).map(aa => {
             ANNSTR.filter{case (ahead,ann) => {ahead ==  aa }}.map{case (ahead,ann) => {ann}}.mkString(",");
          })) else None;
          
          //ANN match {
          //    case Some(annVector) => {
          //      warning("ANN VECTOR: \n      \""+annVector.mkString("\"\n      \"")+"\"","ANN_FIELD_TEST",100);
          //    }
          //    case None => {
          //      //do nothing
          //    }
          //}
          var isBad = false;
          
          val out = altAlleles.zipWithIndex.map{case (alt,altIdx) => {
            var vb = new htsjdk.variant.variantcontext.VariantContextBuilder();
            val altAlleIdx = altIdx + 1;
            vb.loc(vc.getContig(),vc.getStart(),vc.getEnd());
            vb.filters(vc.getFiltersMaybeNull());
            vb.id(vc.getID());
            vb.log10PError(vc.getLog10PError());
            
            val alleleSet = if(alt != Allele.SPAN_DEL) List(refAlle,alt,Allele.SPAN_DEL) else List(refAlle,alt);
            vb.alleles(alleleSet.asJava);
            vb.attribute(vcfCodes.splitAlle_TAG, altAlleles.map(a => {a.getBaseString()}).mkString(",") )
            vb.attribute(vcfCodes.numSplit_TAG, altAlleles.length.toString );
            vb.attribute(vcfCodes.splitIdx_TAG, altIdx.toString );
            
            if(vc.hasGenotypes()){
            vb.genotypes(GenotypesContext.create(
                new java.util.ArrayList[Genotype](vc.getGenotypesOrderedBy(sampNames).asScala.map{gt => {
                  var gb = new GenotypeBuilder(gt.getSampleName());
                  if(gt.hasAD()){
                    val ad = gt.getAD();
                    val refAD = ad.head;
                    val altAD = ad(altAlleIdx);
                    val othAD = ad.zipWithIndex.filter{case (a,i) => {i != 0 && i != altAlleIdx}}.map(_._1).sum
                    gb = gb.AD(Array(refAD,altAD,othAD));
                  }
                  if(gt.hasPL()){
                    val genoPL  = gt.getLikelihoods().getAsPLs();
                    val genoL   = genoPL.map(x => math.pow(10,- x.toDouble / 10.toDouble));
                    gb.PL(Array[Int]( 
                        genoPL(0),
                        genoPL(GenotypeLikelihoods.calculatePLindex(0,altAlleIdx)),
                        genoPL(GenotypeLikelihoods.calculatePLindex(altAlleIdx,altAlleIdx)),
                        (- 10.toDouble * math.log10( multiAllelicIndexStream(numAlle).filter(kp => kp._1 == 0 && kp._2 != altAlleIdx).map(x => genoL(GenotypeLikelihoods.calculatePLindex(x._1,x._2))).sum )).round.toInt,
                        (- 10.toDouble * math.log10( multiAllelicIndexStream(numAlle).filter(kp => (kp._1 == altAlleIdx && kp._2 != altAlleIdx && kp._2 != 0) || (kp._1 != altAlleIdx && kp._2 == altAlleIdx && kp._1 != 0)).map(x => genoL(GenotypeLikelihoods.calculatePLindex(x._1,x._2))).sum )).round.toInt,
                        (- 10.toDouble * math.log10( multiAllelicIndexStream(numAlle).filter(kp => kp._1 != 0 && kp._2 != 0 && kp._1 != altAlleIdx && kp._2 != altAlleIdx).map(x => genoL(GenotypeLikelihoods.calculatePLindex(x._1,x._2))).sum )).round.toInt
                    ));
                  }
                  
                  if(gt.hasDP()) gb = gb.DP(gt.getDP());
                  if(gt.hasGQ()) gb = gb.GQ(gt.getGQ());
                  gb.alleles(gt.getAlleles.asScala.map(a => {
                    if(a.isNoCall()) a;
                    else if(refAlle == a || alt == a) a;
                    else Allele.SPAN_DEL
                  }).asJava); 
                  for((key,attr) <- gt.getExtendedAttributes().asScala){
                    gb.attribute(key,attr);
                  }
                  
                  gb.make();
                }}.toVector.asJava)
            ))}
            simpleCopyAttr.foreach{case (key,attr) => {
              vb.attribute(key,attr);
            }}
            if(clinVarVariants){
              vb.attribute("CLNPROBLEM","0");
              copyAttrCLN.foreach{case(key,attr) => {
                    vb.attribute("ORIG"+key,attr);
                 }}
              val clnAlleList = vc.getAttributeAsList("CLNALLE").asScala.map(x => string2int(x.toString()))
              val clnIdx = clnAlleList.indexWhere(x => x == altAlleIdx);
              if(clnIdx == -1){
                  copyAttrCLN.foreach{case(key,attr) => {
                    vb.attribute(key,"NA");
                  }}
                  warning("CLINVAR CLNALLE TAG IS MALFORMED: CLNALLE="+vc.getAttributeAsString("CLNALLE","?")+" does not contain altAlleIdx = "+altAlleIdx + " (nAlle="+altAlleles.length.toString+")\n"+
                          "       CLNSIG=\""+vc.getAttributeAsList("CLNSIG").asScala.mkString(",")+"\"",
                          "Malformed_ClinVar_CLNALLE_TAG",25);
                  isBad = true;
                  
                  vb.attribute("CLNPROBLEM","1");
              } else {
                //val attrCLN = copyAttrCLN.map{case (key,attr) => {
                copyAttrCLN.foreach{case(key,attr) => {
                  val attrCells = attr.split(",");
                  if(clnIdx >= attrCells.length){
                      warning("CLINVAR "+key+" TAG IS MALFORMED: \""+attr+"\" does not contain clnIdx = "+clnIdx+" (altAlleIdx="+altAlleIdx+", CLNALLE="+vc.getAttributeAsString("CLNALLE","?")+")","Malformed_ClinVar_TAG",25);
                      vb.attribute("CLNPROBLEM","1");
                      isBad = true;
                  } else {
                    vb.attribute(key,attrCells(clnIdx));
                  }
                  if(clnAlleList.length != attrCells.length){
                      warning("CLINVAR "+key+" TAG IS MALFORMED: \""+attr+"\" does not have the same length as CLNALLE="+clnAlleList.mkString(",")+" (altAlleIdx="+altAlleIdx+", CLNALLE="+vc.getAttributeAsString("CLNALLE","?")+")","Malformed_ClinVar_TAG",25);
                      isBad = true;
                  }
                }}
                //}}
              }

            }
            attrA.foreach{case (key,attrArray) => {
              vb.attribute(key,attrArray(altIdx));
            }}
            attrR.foreach{case (key,attrArray) => {
              vb.attribute(key,attrArray(0) + "," + attrArray(altAlleIdx));
            }}
            ANN match {
              case Some(annVector) => {
                vb.attribute("ANN",annVector(altIdx));
              }
              case None => {
                //do nothing
              }
            }
            
            
            /*copyAttrA.foreach{case (key,attr) => {
              val kx = vc.getAttributeAsString(key,"[]").tail.init.split(",")
              vb.attribute(key,kx(altIdx).trim);
            }}
            copyAttrR.foreach{case (key,attr) => {
              val kx = vc.getAttributeAsString(key,"[]").tail.init.split(",")
              vb.attribute(key,kx(0).trim+","+kx(altAlleIdx).trim);
            }}*/
            
            vb.make();
          }}
          if(isBad){
            badCtCt(altAlleles.length) = badCtCt(altAlleles.length) + 1;
          }
          out;
        }
      }), closeAction = () => {
        if(alleCtCt.size > 0){
        Range(1,alleCtCt.keys.max+1).foreach{k => {
          reportln("altCt("+k+"): "+ alleCtCt(k) + " VCF lines","debug");
        }}}
        if(badCtCt.size > 0){
        Range(1,badCtCt.keys.max+1).foreach{k => {
          reportln("badAltCt("+k+"): "+ badCtCt(k) + " VCF lines","debug");
        }}}
      }),newHeader)
    }
  }
  
  //AA,AB,BB,AC,BC,CC
  //0,60,900
  
  class CmdRecodeClinVarCLN extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "RecodeClinVarCLN", 
          quickSynopsis = "", 
          synopsis = "", 
          description = " " + ALPHA_WARNING,   
          argList = 
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "..."
                                        ) ::                     
                    new FinalArgument[String](
                                         name = "invcf",
                                         valueName = "variants.vcf",
                                         argDesc = "infput VCF file" // description
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
       
       VcfAnnotateTX.RecodeClinVarCLN(
                   chromList = parser.get[Option[List[String]]]("chromList")
                   ).walkVCFFile(
                   infile    = parser.get[String]("invcf"),
                   outfile   = parser.get[String]("outvcf"),
                   chromList = parser.get[Option[List[String]]]("chromList")
                   )
     }   
    }
  }
  
  case class RecodeClinVarCLN(
                 chromList : Option[List[String]],
                 vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ) extends internalUtils.VcfTool.VCFWalker {
    
    reportln("Creating recodeClinVarCLN()","note")

    val infoCLN = Set[String]("CLNHGVS","CLNSRC","CLNORIGIN","CLNSRCID","CLNSIG","CLNDSDB","CLNDSDBID","CLNREVSTAT","CLNACC","CLNDBN");
    val infoCLNALL = infoCLN ++ Set[String]("CLNALLE");
    
    val alleCtCt = new scala.collection.mutable.HashMap[Int, Int]().withDefaultValue(0)
    val badCtCt = new scala.collection.mutable.HashMap[Int, Int]().withDefaultValue(0)
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      val oldHeaderLines =  vcfHeader.getMetaDataInInputOrder().asScala.toList;
      
      val (keepOldLines,oldClnLines) = oldHeaderLines.partition(headerline => {
        headerline match {
          case h : VCFInfoHeaderLine => {
            ! infoCLNALL.contains(h.getID()); 
          }
          case _ => true;
        }
      })
      val oldClnHeader = oldClnLines.map(headerLine => {
        headerLine match {
          case h : VCFInfoHeaderLine => h;
          case _ => null;
        }
      })
       

      val newHeaderLines = oldClnHeader.map(h => {
            new VCFInfoHeaderLine(h.getID(), VCFHeaderLineCount.A, h.getType(), h.getDescription())
          }) ++ List(
            new VCFInfoHeaderLine("CLNERROR",     1, VCFHeaderLineType.String, "")
          ) ++ oldClnHeader.map(h => {
            new VCFInfoHeaderLine("ORIG"+h.getID(),     VCFHeaderLineCount.UNBOUNDED, h.getType(), h.getDescription() + " (ORIGINAL UNFIXED VALUE FROM CLINVAR)")
          })
      
      val replacementHeaderLines = keepOldLines ++ newHeaderLines;
          
      val newHeader = internalUtils.VcfTool.replaceHeaderLines(vcfHeader,replacementHeaderLines);
      
      reportln("Walking input VCF...","note")
      return ( 
      addIteratorCloseAction(vcIter.map(v => {
        runRecodeClinVarCLN(v);
      }), closeAction = () => {
        if(alleCtCt.size > 0){
        Range(1,alleCtCt.keys.max+1).foreach{k => {
          reportln("altCt("+k+"): "+ alleCtCt(k) + " VCF lines","debug");
        }}}
        if(badCtCt.size > 0){
        Range(1,badCtCt.keys.max+1).foreach{k => {
          reportln("badAltCt("+k+"): "+ badCtCt(k) + " VCF lines","debug");
        }}}
      }), newHeader );
    }
    
    def runRecodeClinVarCLN(v : VariantContext) : VariantContext = {
      var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
      vb.attribute("CLNERROR","OK");
      var isErr = false;
      val clnAlle = v.getAttributeAsList("CLNALLE").asScala.map(x => string2int(x.toString())).toVector;
          val alleles = internalUtils.VcfTool.getAllelesInOrder(v).toVector;
          val numAlle = alleles.length;
          val refAlle = alleles.head;
          val altAlleles = alleles.tail;
          alleCtCt(altAlleles.length) = alleCtCt(altAlleles.length) + 1;
      val clnMap = Range(1,alleles.length).map(alleIdx => {
        clnAlle.indexWhere(clnAlleNum => {
          clnAlleNum == alleIdx;
        })
      }).toVector
      
      
      
      infoCLN.foreach(key => {
        val attr = v.getAttributeAsList(key).asScala.map(_.toString());
        vb.attribute("ORIG"+key, attr.mkString(","));
        if(attr.length != clnAlle.length){
          isErr = true;
          vb.attribute("CLNERROR","ERR");
        } else {
          val vcfAllelesNotFound = clnMap.count(_ == -1);
          val clnAllelesNotFound = clnAlle.count(i => { i == -1 });
          
          if(vcfAllelesNotFound > 0 && clnAllelesNotFound > 0){
              isErr = true;
              vb.attribute("CLNERROR","ERR");
          }
          
          val newAttr = clnMap.map(clnIdx => {
            if(clnIdx == -1){
              //isErr = true;
              //vb.attribute("CLNERROR","ERR");
              "NA";
            } else {
              attr(clnIdx);
            }
          });
          vb.attribute(key, newAttr.mkString(","));
        }
      })
      vb.attribute("ORIGCLNALLE",v.getAttributeAsList("CLNALLE").asScala.map(_.toString()).mkString(","));
      vb.attribute("CLNALLE",Range(0,altAlleles.length).map(x => (x + 1).toString()).mkString(","));
      if(isErr){
        badCtCt(altAlleles.length) = badCtCt(altAlleles.length) + 1;
      }
      
      return vb.make();
    }
    
    reportln("recodeClinVarCLN() Created...","note")
    
  }
  

  
  
}




























