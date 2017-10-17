package internalTests

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

import htsjdk.samtools._;
import internalUtils.fileUtils._;
import java.io.File._;
import scala.util.Random._;

import internalUtils.genomicAnnoUtils.GenomicArrayOfSets
import scala.collection.JavaConverters._

object makeBedWiggle {
   
  class CmdCodingCoverageStats extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "CodingCoverageStats", 
          quickSynopsis = "", 
          synopsis = "", 
          description = ""+ALPHA_WARNING,
          argList = 
                    new UnaryArgument( name = "infileList",
                                         arg = List("--infileList"), // name of value
                                         argDesc = "If this option is used, then instead of a single input file the input file(s) will be assumed "+
                                                   "to be a file containing a list of input files to parse in order. If multiple VCF files are specified, "+
                                                   "the vcf lines will be concatenated and the header will be taken from the first file."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument(   name = "singleEnded", 
                                         arg = List("--singleEnded"), // name of value
                                         argDesc = "Flag for single-end data. Note that many other options do not apply in this case (for example: option --countPairsTogether does nothing in single-end mode)" 
                                       ) ::
                    new BinaryOptionArgument[String](
                                         name = "onTargetBed", 
                                         arg = List("--onTargetBed"), 
                                         valueName = "targetregion.bed",  
                                         argDesc =  ""
                                        ) ::
                    new BinaryArgument[String](
                                         name = "trackTitle", 
                                         arg = List("--trackTitle"), 
                                         valueName = "bpCoverage",  
                                         argDesc =  "",
                                         defaultValue = Some("bpCoverage")
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile.bam",
                                         argDesc = "Input bam file."// description
                                        ) ::
                    new FinalArgument[String](
                                         name = "inputSavedTxFile",
                                         valueName = "txdata.txt.gz",
                                         argDesc = "Input txdata file."// description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfileprefix",
                                         valueName = "outfile",
                                         argDesc = "The output file prefix. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         CodingCoverageStats(
             infile = parser.get[String]("infile"),
             outfileprefix = parser.get[String]("outfileprefix"),
             isSingleEnd = parser.get[Boolean]("singleEnded"),
             onTargetBed = parser.get[Option[String]]("onTargetBed"),
             inputSavedTxFile = parser.get[String]("inputSavedTxFile"),
             trackTitle = parser.get[String]("trackTitle")
         ).run();
       }
     }
  }
  
  case class CodingCoverageStats(infile : String, outfileprefix :String, isSingleEnd : Boolean, 
                                 onTargetBed : Option[String], inputSavedTxFile : String, 
                                 trackTitle : String = "bpcoverage",
                                 windowCt : Int = 2000,
                                 spannedWindowSize : Int = 1000,
                                 chromLengthFile : Option[String] = None){
    
    val BED_FILE_INTERNAL_TAG = "ON_TARGET_BED_FILE"
    
    val coverageThresholds = Vector[(Int,Int)]((0,1),(1,2),(2,3),(3,4),(4,5),(5,10),(10,20),(20,30),(30,40),(40,50),(50,60),(60,70),(70,80),(80,90),(90,100),(100,Integer.MAX_VALUE));
    //var coverageCounts : Vector[(String,Array[Int])] = Vector[(String,Array[Int])]()
    
    val onTargetFilter = onTargetBed match {
      case Some(bedfile) => {
        (txset : Set[String]) => txset.contains(BED_FILE_INTERNAL_TAG) && txset.size >= 2;
      }
      case None => {
        (txset : Set[String]) => txset.size >= 1;
      }
    }
    val targetBedIv = onTargetBed match {
      case Some(bedfile) => {
        Some(getLinesSmartUnzip(bedfile).map{line => {
          val cells = line.split("\t");
          (cells(0),string2int(cells(1)),string2int(cells(2)))
        }}.toVector)
      }
      case None => {
        None;
      }
    }
    val chromLens : Option[Map[String,Int]] = chromLengthFile match {
      case Some(clf) => {
        Some(getLinesSmartUnzip(clf).map{line => {
          val cells = line.split("\t");
          (cells(0),string2int(cells(1)))
        }}.toMap)
      }
      case None => None;
    }
    
    val (iterRaw, header, attrib) : (Iterator[(SAMRecord,SAMRecord)],htsjdk.samtools.SAMFileHeader,internalUtils.SamTool.SamFileAttributes) = 
      internalUtils.SamTool.getPairedEndReader(inbam = infile, 
                         isSingleEnd = isSingleEnd,
                         sortPairsByFirstPosition = true,
                         unsorted = true,
                         maxPhredScore = 41,
                         stopAfterNReads = None,
                         testRunLineCt  = 10000,
                         testRun = false,
                         maxReadLength = None)
    val iter = iterRaw.buffered;
    
    val TXSeq : Vector[internalUtils.TXUtil] = {
          val ipr = internalUtils.stdUtils.IteratorProgressReporter_ThreeLevel("lines", 200, 1000 , 2000 )
          val wrappedIter = internalUtils.stdUtils.wrapIteratorWithProgressReporter(getLinesSmartUnzip(inputSavedTxFile) , ipr )
          //val txs = scala.collection.mutable.AnyRefMap[String,internalUtils.TXUtil]();
          wrappedIter.map{line => {
            val tx = internalUtils.TXUtil.buildTXUtilFromString(line);
            tx
          }}.toVector
    }
    
    val totalCounts_tx = Array.fill(coverageThresholds.length)(0);
    val totalCounts_cd = Array.fill(coverageThresholds.length)(0);  
    
    case class ChromDataHolder(chr : String){
      var txArray = GenomicArrayOfSets[String](false);
      var codingArray = GenomicArrayOfSets[String](false);
      TXSeq.withFilter(tx => {tx.chrom == chr} ).foreach{ tx => {
        tx.gSpans.foreach{ case (start,end) => {
                txArray.addSpan(GenomicInterval(chromName = tx.chrom, strand = '.', start = start, end = end), tx.txID);
        }}
        tx.gSpansCDS.foreach{ case (start,end) => {
                codingArray.addSpan(GenomicInterval(chromName = tx.chrom, strand = '.', start = start, end = end), tx.txID);
        }}
      }}
      targetBedIv match {
        case Some(tbiv) => {
          tbiv.withFilter{ case (chrom,start,end) => { chrom == chr}}.foreach{ case (chrom,start,end) => {
            val iv = GenomicInterval(chromName = chrom, strand = '.', start = start, end = end)
            txArray.addSpan(iv, BED_FILE_INTERNAL_TAG);
            codingArray.addSpan(iv,BED_FILE_INTERNAL_TAG);
          }}
        }
        case None => {
          //do nothing!
        }
      }
      txArray.finalizeStepVectors;
      codingArray.finalizeStepVectors;
      var stepCountArrays = txArray.getSteps(chr,'.').withFilter{ case (iv,txset) => onTargetFilter(txset)}.map{ case (iv,txset) => {
        (iv,(txset,Array.fill(iv.end - iv.start)(0)))
      }}.toMap;
      var codingStepCountArrays = codingArray.getSteps(chr,'.').withFilter{ case (iv,txset) => onTargetFilter(txset)}.map{ case (iv,txset) => {
        (iv,(txset,Array.fill(iv.end - iv.start)(0)))
      }}.toMap;
      var ivlist_tx =  txArray.getSteps(chr,'.').withFilter{ case (iv,txset) => onTargetFilter(txset)}.map{ case (iv,txset) => iv }.toList
      var ivlist_cd =  codingArray.getSteps(chr,'.').withFilter{ case (iv,txset) => onTargetFilter(txset)}.map{ case (iv,txset) => iv }.toList
      var chromWindowSums_cds = Array.fill(windowCt)(0);
      var chromSize = ivlist_cd.map(iv => iv.end - iv.start).sum;
      var windowSize = chromSize / windowCt + 1;
      
      var spannedWindowSums_cds = Array.fill(chromSize / spannedWindowSize + 1)(0)
      /*var chromWindowSums_cds = ivlist_cd.tail.foldLeft(Vector(ivlist_cd.head)){case (soFar,iv) => {
        if(soFar.last.end == iv.start){
          soFar.init :+ GenomicInterval(chromName = iv.chromName, strand = '.', start = soFar.last.start, end = iv.end)
        } else {
          soFar :+ iv
        }
      }}.map{ iv => {
        
      }}*/
      
      def writeChrom(txOut : internalUtils.fileUtils.WriterUtil,
                     cdOut : internalUtils.fileUtils.WriterUtil,
                     summaryOutTx : internalUtils.fileUtils.WriterUtil,
                     summaryOutCd : internalUtils.fileUtils.WriterUtil,
                     cdWindows : internalUtils.fileUtils.WriterUtil,
                     spannedCdWindows :  internalUtils.fileUtils.WriterUtil){
            reportln("writing chrom: " + chr + " [" + getDateAndTimeString+"]","note");
            val currChromCounts_tx = Array.fill(coverageThresholds.length)(0);
            val currChromCounts_cd = Array.fill(coverageThresholds.length)(0);
            
            ivlist_tx.map{ iv => (iv, stepCountArrays(iv)) }.foreach{ case (iv,(txset,countArray)) => {
              txOut.write("#"+"\t"+iv.chromName+"\t"+iv.start +"\t"+iv.end+"\t"+txset.filter(tx => tx != BED_FILE_INTERNAL_TAG).toVector.sorted.mkString(",")+"\n");
              txOut.write("fixedStep chrom="+iv.chromName+" start="+(iv.start+1)+" step=1\n");
              txOut.write(countArray.mkString("\n")+"\n");
              countArray.foreach{ct => {
                coverageThresholds.zipWithIndex.foreach{ case ((ts,te),i) => { if(ct < te && ct >= ts) currChromCounts_tx(i) += 1 }}
              }}
            }}
            var currPos = 0;
            ivlist_cd.map{ iv => (iv, codingStepCountArrays(iv)) }.foreach{ case (iv,(txset,countArray)) => {
              cdOut.write("#"+"\t"+iv.chromName+"\t"+iv.start +"\t"+iv.end+"\t"+txset.filter(tx => tx != BED_FILE_INTERNAL_TAG).toVector.sorted.mkString(",")+"\n");
              cdOut.write("fixedStep chrom="+iv.chromName+" start="+(iv.start+1)+" step=1\n");
              cdOut.write(countArray.mkString("\n")+"\n");
              countArray.foreach{ct => {
                coverageThresholds.zipWithIndex.foreach{ case ((ts,te),i) => { if(ct < te && ct >= ts) currChromCounts_cd(i) += 1 }}
                chromWindowSums_cds(currPos / windowSize) += ct;
                spannedWindowSums_cds(currPos / spannedWindowSize) += ct;
                currPos += 1;
              }}
            }}
            
            currChromCounts_tx.indices.foreach{i => totalCounts_tx(i) += currChromCounts_tx(i)}
            currChromCounts_cd.indices.foreach{i => totalCounts_cd(i) += currChromCounts_cd(i)}
            summaryOutTx.write(chr+"\t"+currChromCounts_tx.mkString("\t")+"\n");
            summaryOutCd.write(chr+"\t"+currChromCounts_cd.mkString("\t")+"\n");
            cdWindows.write(chr+"\t"+chromWindowSums_cds.map{ct => ct.toDouble / windowSize.toDouble}.mkString("\t")+"\n");
            spannedWindowSums_cds.zipWithIndex.foreach{ case (ct,idx) => {
              spannedCdWindows.write(chr+"\t"+idx+"\t"+ct.toDouble / spannedWindowSize.toDouble+"\n")
            }}
            //spannedWindowSums_cds
            summaryOutTx.flush();
            summaryOutCd.flush();
            cdWindows.flush();
            spannedCdWindows.flush();
      }
    }
    var currChrom : String = iter.head._1.getContig();
    var cdata = ChromDataHolder(currChrom);

    def run(){
      reportln("Initializing output files ["+getDateAndTimeString+"]","note");
      val txOut = openWriterSmart(outfileprefix + "gene.baseDepths.txt.gz");
      val cdOut = openWriterSmart(outfileprefix + "cds.baseDepths.txt.gz");
      txOut.write("track type=wiggle_0 name="+trackTitle+"_Genic\n");
      cdOut.write("track type=wiggle_0 name="+trackTitle+"_CDS\n");
      val summaryOutTx = openWriterSmart(outfileprefix + "gene.depthSummaryByChrom.txt");
      val summaryOutCd = openWriterSmart(outfileprefix + "cds.depthSummaryByChrom.txt");
      summaryOutTx.write("chrom"+"\t"+coverageThresholds.init.map{ case (ts,te) => if(ts+1==te) ts else ts+"to"+(te-1)}.mkString("\t")+"\tge"+coverageThresholds.last._1+"\n");
      summaryOutCd.write("chrom"+"\t"+coverageThresholds.init.map{ case (ts,te) => if(ts+1==te) ts else ts+"to"+(te-1)}.mkString("\t")+"\tge"+coverageThresholds.last._1+"\n");
      
      val cdWindows = openWriterSmart(outfileprefix + "cds.windowedDepths.equalNumWindows.txt");
      cdWindows.write("chrom\t"+Range(0,windowCt).mkString("\t")+"\n");
      val spannedCdWindows = openWriterSmart(outfileprefix + "cds.windowedDepths.equalSizeWindows.txt");
      spannedCdWindows.write("chrom\tindex\tct\n");
      
      reportln("Starting iteration ["+getDateAndTimeString+"]","note");
      for((r1,r2) <- iter){
        if((! r1.getReadUnmappedFlag()) && (! r1.getReadFailsVendorQualityCheckFlag()) && (! r2.getReadUnmappedFlag()) && (! r2.getReadFailsVendorQualityCheckFlag()) && r1.getAlignmentBlocks().size() != 0 && r2.getAlignmentBlocks().size() != 0){
          if(r1.getContig() != currChrom){
            /////////////////////////////////////////////////////////////////////////
            cdata.writeChrom(txOut=txOut ,  cdOut=cdOut , summaryOutTx=summaryOutTx, summaryOutCd=summaryOutCd ,cdWindows =cdWindows,spannedCdWindows=spannedCdWindows)
            /////////////////////////////////////////////////////////////////////////
            reportln("Switching to chrom: " + r1.getContig()+ " [" + getDateAndTimeString+"]","note");
            currChrom = r1.getContig();
            cdata = ChromDataHolder(currChrom);
          }
          
          val blocks = getOverlappedPairBlocks(r1,r2);
          blocks.foreach{ case (start,end) => {
            val iv = GenomicInterval(chromName = currChrom,strand='.',start=start,end=end);
            val txSteps = cdata.txArray.findIntersectingSteps(iv);
            val cdSteps = cdata.codingArray.findIntersectingSteps(iv);
            txSteps.withFilter{ case (iv,txset) => onTargetFilter(txset)}.foreach{ case (stepIV,stepTxSet) => {
              val stepCts = cdata.stepCountArrays(stepIV)._2;
              val from = if(iv.start > stepIV.start) iv.start - stepIV.start else 0;
              val to = if(iv.end < stepIV.end) stepCts.length - (stepIV.end - iv.end) else stepCts.length;
              Range(from,to).foreach{ i => {
                  stepCts(i) += 1;
              }}
            }}
            cdSteps.withFilter{ case (iv,txset) => onTargetFilter(txset)}.foreach{ case (stepIV,stepTxSet) => {
              val stepCts = cdata.codingStepCountArrays(stepIV)._2;
              val from = if(iv.start > stepIV.start) iv.start - stepIV.start else 0;
              val to = if(iv.end < stepIV.end) stepCts.length - (stepIV.end - iv.end) else stepCts.length;
              Range(from,to).foreach{ i => {
                  stepCts(i) += 1;
              }}
            }}
          }}
        }
      }
      cdata.writeChrom(txOut=txOut ,  cdOut=cdOut , summaryOutTx=summaryOutTx, summaryOutCd=summaryOutCd ,cdWindows =cdWindows,spannedCdWindows=spannedCdWindows)
      
      
      reportln("Finished iteration ["+getDateAndTimeString+"]","note");

      summaryOutTx.write("TOTAL"+"\t"+totalCounts_tx.mkString("\t")+"\n");
      summaryOutTx.write("PCT"+"\t"+totalCounts_tx.map{cts => cts.toDouble / totalCounts_tx.sum.toDouble}.mkString("\t")+"\n");
      val cumsumtx = totalCounts_tx.reverse.scanLeft(0)( (soFar,curr) => soFar + curr ).tail.reverse
      summaryOutTx.write("Cumulative"+"\t"+cumsumtx.mkString("\t")+"\n");
      summaryOutTx.write("CumulativePct"+"\t"+cumsumtx.map{_.toDouble / totalCounts_tx.sum.toDouble}.mkString("\t")+"\n");
      
      summaryOutCd.write("TOTAL"+"\t"+totalCounts_cd.mkString("\t")+"\n");
      summaryOutCd.write("PCT"+"\t"+totalCounts_cd.map{cts => cts.toDouble / totalCounts_cd.sum.toDouble}.mkString("\t")+"\n");
      val cumsumcd = totalCounts_cd.reverse.scanLeft(0)( (soFar,curr) => soFar + curr ).tail.reverse
      summaryOutTx.write("Cumulative"+"\t"+cumsumcd.mkString("\t")+"\n");
      summaryOutTx.write("CumulativePct"+"\t"+cumsumcd.map{_.toDouble / totalCounts_cd.sum.toDouble}.mkString("\t")+"\n");
      
      cdWindows.close(); 
      summaryOutTx.close();
      summaryOutCd.close();
      txOut.close();
      cdOut.close();
      spannedCdWindows.close();
    } //end run() method
  }
  
  def getReadBlocks(r : SAMRecord) : Vector[(Int,Int)] = {
    r.getAlignmentBlocks().asScala.toVector.map((block) => {
      (block.getReferenceStart() - 1, block.getReferenceStart() - 1 + block.getLength());
    });
  }
  
  def getOverlappedPairBlocks(r1 : SAMRecord, r2 : SAMRecord) : Vector[(Int,Int)] = {
    val r1b = getReadBlocks(r1);
    val r2b = getReadBlocks(r2);
    
    //def blocksOverlap(b1 : (Int,Int), b2 : (Int,Int)) : Boolean = {
    //  b1._1 <= b2._2 && b2._1 <= b1._2;
    //}
    //def mergeBlocks(b1 : (Int,Int), b2 : (Int,Int)) : (Int,Int) = {
    //  (math.min(b1._1,b2._1), math.max(b1._2, b2._2))
    //}
    //val r2bOverlap : Seq[(Int,Int)] = r2b.filter((b2) => r1b.exists(blocksOverlap(_,b2)));
    //val r2bNonOverlap : Seq[(Int,Int)] = r2b.filterNot((b2) => r1b.exists(blocksOverlap(_,b2)));
    val merged = (r1b ++ r2b).sorted
    merged.tail.foldLeft(Vector(merged.head))( (soFar,curr) =>{
      if(curr._1 <= soFar.last._2){
        soFar.updated(soFar.length - 1, (soFar.last._1, math.max(curr._2, soFar.last._2)));
      } else {
        soFar :+ curr;
      }
    })
  }
}


























