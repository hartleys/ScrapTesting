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
  
  class compareVcfs extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "compareVcfs", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryArgument[String](
                                         name = "GenoTag1", 
                                         arg = List("--GenoTag1"), 
                                         valueName = "GT",  
                                         argDesc =  "",
                                         defaultValue = Some("GT")
                                        ) ::
                    new BinaryArgument[String](
                                         name = "GenoTag2", 
                                         arg = List("--GenoTag2"), 
                                         valueName = "GT",  
                                         argDesc =  "",
                                         defaultValue = Some("GT")
                                        ) ::
                    new UnaryArgument( name = "infileList",
                                         arg = List("--infileList"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "noGzipOutput",
                                         arg = List("--noGzipOutput"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::
                    new FinalArgument[String](
                                         name = "infile1",
                                         valueName = "variants1.vcf",
                                         argDesc = "input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile2",
                                         valueName = "variants2.vcf",
                                         argDesc = "input VCF file. Can be gzipped or in plaintext." // description
                                        ) :: 
                    new FinalArgument[String](
                                         name = "outfileprefix",
                                         valueName = "outfileprefix",
                                         argDesc = "The output file."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         VcfStreamCompare(
                          infile1 = parser.get[String]("infile1"),
                          infile2 = parser.get[String]("infile2"),
                          outfile = parser.get[String]("outfileprefix"),
                          chromList = parser.get[Option[List[String]]]("chromList"),
                          genoTag1 = parser.get[String]("GenoTag1"),
                          genoTag2 = parser.get[String]("GenoTag2"),
                          infileList = parser.get[Boolean]("infileList"),
                          gzipOutput = ! parser.get[Boolean]("noGzipOutput")
         ).run()
       }
     }
    
  }
  

  
  
  case class VcfStreamCompare(infile1 : String, infile2 : String, outfile : String,
                              chromList : Option[List[String]],
                              genoTag1 : String,
                              genoTag2 : String,
                              file2tag : String = "B_",
                              file2Desc : String = "(For Alt Build) ",
                              infileList : Boolean,
                              gzipOutput : Boolean){
    val (vcIterRaw1, vcfHeader1) = getSVcfIterators(infile1,chromList,None,inputFileList = infileList);
    val (vcIterRaw2, vcfHeader2) = getSVcfIterators(infile2,chromList,None,inputFileList = infileList, withProgress = false);
    val (vcIter1,vcIter2) = (vcIterRaw1.buffered, vcIterRaw2.buffered)
    var currChrom = vcIter1.head.chrom;
    
    val outfileSuffix = if(gzipOutput) ".gz" else "";
    
    val matchIdx = vcfHeader1.titleLine.sampleList.zipWithIndex.flatMap{ case (s,idx1) => {
      val idx2 = vcfHeader2.titleLine.sampleList.indexOf(s)
      if(idx2 == -1){
        None
      } else {
        Some((s,idx1,idx2))
      }
    }}
    val matchIdxIdx = matchIdx.indices;
    /*
     * A>C, T>G
     * A>T, T>A
     * A>G, T>C
     * C>A, G>T
     * C>T, G>A
     * C>G, G>C
     */
    val baseSwapList = Seq( (("A","C"),("T","G")),
                            (("A","T"),("T","A")),
                            (("A","G"),("T","C")),
                            (("C","A"),("G","T")),
                            (("C","T"),("G","A")),
                            (("C","G"),("G","C"))
                          );
    val fGenoTag2 = file2tag + genoTag2;

    val vcfHeaderOut = SVcfHeader(infoLines = vcfHeader1.infoLines, 
                        formatLines = vcfHeader1.formatLines ++ vcfHeader2.formatLines.map{ fl => {
                          SVcfCompoundHeaderLine(in_tag = fl.in_tag, ID = file2tag+fl.ID, Number = fl.Number, Type = fl.Type, desc = file2Desc + fl.desc)
                        }},
                        otherHeaderLines = vcfHeader1.otherHeaderLines,
                        titleLine = SVcfTitleLine(sampleList = matchIdx.map{_._1}.toSeq)
                     )

    
    def isFirst : Boolean = vcIter1.hasNext && (
                        (! vcIter2.hasNext) || 
                        ( vcIter1.head.chrom == vcIter2.head.chrom && vcIter1.head.pos <= vcIter2.head.pos) ||
                        ( vcIter1.head.chrom == currChrom && vcIter2.head.chrom != currChrom ) 
                   );
    
    val outM = openWriterSmart(outfile+"shared.vcf"+outfileSuffix);
    //val outMM = openWriterSmart(outfile+"sharedMis.vcf.gz");
    val out1 = openWriterSmart(outfile+"F1.vcf"+outfileSuffix);
    val out2 = openWriterSmart(outfile+"F2.vcf"+outfileSuffix);
    val writers = Seq(outM,out1,out2);
    
    val matchWriter = openWriterSmart(outfile+"summary.match.txt"+outfileSuffix);
    val AWriter = openWriterSmart(outfile+"summary.A.txt"+outfileSuffix);
    val BWriter = openWriterSmart(outfile+"summary.B.txt"+outfileSuffix);

    writers.foreach(out => {
      vcfHeaderOut.getVcfLines.foreach{ l =>
        out.write(l+"\n");
      }
    })
    
    var matchCt = 0;
    var at1not2ct = 0;
    var at2not1ct = 0;
    var alleAt1not2ct = 0;
    var alleAt2not1ct = 0;
    
    def run(){
      reportln("> ITERATE()","debug");
      while(vcIter1.hasNext || vcIter2.hasNext){
        reportln("> ITERATE()","debug");
        iterate();
      }
      writers.foreach(out => {
        out.close();
      })
      matchWriter.write(
            "sample.ID\t"+matchCountFunctionList_FINAL.map{ case (id,arr,varFcn,fcn) => {
              id
            }}.mkString("\t")+"\n"
      );
      AWriter.write(
            "sample.ID\t"+mmCountFunctionList_A.map{ case (id,arr,varFcn,fcn) => {
              id
            }}.mkString("\t")+"\n"
      );
      BWriter.write(
            "sample.ID\t"+mmCountFunctionList_B.map{ case (id,arr,varFcn,fcn) => {
              id
            }}.mkString("\t")+"\n"
      );
      
      
      matchIdx.zipWithIndex.foreach{ case ((sampid,idx1,idx2),i) => {
        matchWriter.write(
            sampid+"\t"+matchCountFunctionList_FINAL.map{ case (id,arr,varFcn,fcn) => {
              arr(i);
            }}.mkString("\t")+"\n"
        );
        AWriter.write(
            sampid+"\t"+mmCountFunctionList_A.map{ case (id,arr,varFcn,fcn) => {
              arr(i);
            }}.mkString("\t")+"\n"
        );
        BWriter.write(
            sampid+"\t"+mmCountFunctionList_B.map{ case (id,arr,varFcn,fcn) => {
              arr(i);
            }}.mkString("\t")+"\n"
        );
      }}
      matchWriter.close();
      AWriter.close();
      BWriter.close();
      
    }
    
    def iterate(){
      val vcos = if(
                    (vcIter1.hasNext & ! vcIter2.hasNext)  || (vcIter2.hasNext & ! vcIter1.hasNext) || 
                    (vcIter1.head.pos != vcIter2.head.pos) || (vcIter1.head.chrom != vcIter2.head.chrom)){
        reportln("   Iterating Isolated position.","deepDebug");
        if(isFirst){
          val vc = vcIter1.next();
          reportln("   Iterating VC1 Isolated Position: "+vc.chrom+":"+vc.pos,"debug");
          at1not2ct += 1;
          writeVC1(vc);
        } else {
          val vc = vcIter2.next().getOutputLine();
          reportln("   Iterating VC2 Isolated Position: "+vc.chrom+":"+vc.pos,"debug");
          at2not1ct += 1;
          writeVC2(vc);
        };
      } else {
        reportln("   Iterating Shared Position: "+vcIter1.head.chrom+":"+vcIter1.head.pos,"debug");
        iteratePos();
      }
    }
    
    def isCalled(vc : SVcfVariantLine, idx : Int, i : Int) : Boolean = {
      ! vc.genotypes.genotypeValues(idx)(i).contains('.')
    }
    def isCalled(vc : SVcfVariantLine, idx1 : Int, idx2 : Int, i : Int) : Boolean = {
      isCalled(vc,idx1,i) && isCalled(vc,idx2,i);
    }
    def isMatch(vc : SVcfVariantLine,idx1 : Int, idx2: Int, i : Int) : Boolean = {
      vc.genotypes.genotypeValues(idx1)(i) == vc.genotypes.genotypeValues(idx2)(i) && isCalled(vc,idx1,i)
    }
    def isMisMatch(vc : SVcfVariantLine,idx1 : Int, idx2: Int, i : Int) : Boolean = {
      vc.genotypes.genotypeValues(idx1)(i) != vc.genotypes.genotypeValues(idx2)(i) && isCalled(vc,idx1,i) && isCalled(vc,idx2,i)
    }
    def isSNV(vc : SVcfVariantLine) : Boolean = {
      vc.ref.length == vc.alt.length && vc.ref.length == 1
    }
    def isRef(vc : SVcfVariantLine, idx : Int, i : Int) : Boolean = {
       vc.genotypes.genotypeValues(idx)(i) == "0/0";
    }
    def isHomAlt(vc : SVcfVariantLine, idx : Int, i : Int) : Boolean = {
      vc.genotypes.genotypeValues(idx)(i) == "1/1";
    }
    def isHet(vc : SVcfVariantLine, idx : Int, i : Int) : Boolean = {
      vc.genotypes.genotypeValues(idx)(i) == "0/1";
    }
    def isAnyAlt(vc : SVcfVariantLine, idx : Int, i : Int) : Boolean = {
      isHet(vc,idx,i) || isHomAlt(vc,idx,i);
    }
    def isOther(vc : SVcfVariantLine, idx : Int, i : Int) : Boolean = {
      vc.genotypes.genotypeValues(idx)(i).length != 3 || vc.genotypes.genotypeValues(idx)(i).charAt(0) != '2' || vc.genotypes.genotypeValues(idx)(i).charAt(2) == '2'
    }
    
    val matchCountFunctionList_BASE : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int,Int) => Boolean)] = Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int,Int) => Boolean)](
        ("noCall",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              ! isCalled(vc,idx1,i) && ! isCalled(vc,idx2,i)
        }),
        ("calledA",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isCalled(vc,idx1,i) && ! isCalled(vc,idx2,i)
        }),
        ("calledB",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              ! isCalled(vc,idx1,i) && isCalled(vc,idx2,i)
        }),
        ("called",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isCalled(vc,idx1,i) && isCalled(vc,idx2,i)
        }),
        ("match",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMatch(vc,idx1,idx2,i)
        }),
        ("mismatch",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i)
        }),
        ("mismatch.Ref.Alt",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isRef(vc,idx1,i) && isAnyAlt(vc,idx2,i);
        }),
        ("mismatch.Alt.Ref",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isRef(vc,idx2,i) && isAnyAlt(vc,idx1,i);
        }),
        ("mismatch.Alt.Alt",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && (! isRef(vc,idx1,i)) && (! isRef(vc,idx2,i));
        }),
        ("mismatch.Ha.Het",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isHomAlt(vc,idx1,i) && isHet(vc,idx2,i);
        }),
        ("mismatch.Het.Ha",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isHomAlt(vc,idx2,i) && isHet(vc,idx1,i);
        }),
        ("mismatch.Ref.Het",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isRef(vc,idx1,i) && isHet(vc,idx2,i);
        }),
        ("mismatch.Het.Ref",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isRef(vc,idx2,i) && isHet(vc,idx1,i);
        }),
        ("mismatch.Ref.Ha",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isRef(vc,idx1,i) && isHomAlt(vc,idx2,i);
        }),
        ("mismatch.Ha.Ref",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && isRef(vc,idx2,i) && isHomAlt(vc,idx1,i);
        }),
        ("mismatch.Other",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx1: Int, idx2 :Int, i : Int) => {
              isMisMatch(vc,idx1,idx2,i) && (isOther(vc,idx1,i) || isOther(vc,idx2,i))
        })
    );
    val matchCountFunctionList_BYTYPE = matchCountFunctionList_BASE.map{ case (id,arr,varFcn,fcn) => {
      (id + "_SNV",Array.fill[Int](matchIdx.length)(0),(vc : SVcfVariantLine) => isSNV(vc),fcn)
    }} ++ matchCountFunctionList_BASE.map{ case (id,arr,varFcn,fcn) => {
      (id + "_INDEL",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => ! isSNV(vc),fcn)
    }}
    
    val matchCountFunctionList_BYSWAP = matchCountFunctionList_BASE.flatMap{ case (id,arr,varFcn,fcn) => {
      baseSwapList.map{ case ((r1,a1),(r2,a2)) => {
        (id+"_SNV_SWAP."+r1+a1,Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => {
          isSNV(vc) && ((vc.ref == r1 && vc.alt == a1) || (vc.ref == r2 && vc.alt == a2))
        },fcn)
      }}
    }}
    
    val matchCountFunctionList_BASE2 : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int,Int) => Boolean)] = matchCountFunctionList_BASE ++ matchCountFunctionList_BYTYPE ++ matchCountFunctionList_BYSWAP;
    val matchCountFunctionList_BASE2_CFILTSET = matchCountFunctionList_BASE2.map{ case (id,arr, varFcn,fcn) => {
      ("CPASSAB_"+id,Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => {
        varFcn(vc) && vc.filter == "."
      },fcn)
    }} ++ matchCountFunctionList_BASE2.map{ case (id,arr, varFcn,fcn) => {
      ("CPASSA_"+id,Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => {
        varFcn(vc) && vc.filter != "." && vc.filter.split(",")(0) == ".";
      },fcn)
    }} ++ matchCountFunctionList_BASE2.map{ case (id,arr, varFcn,fcn) => {
      ("CPASSB_"+id,Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => {
        varFcn(vc) && vc.filter != "." && vc.filter.split(",")(1) == ".";
      },fcn)
    }} ++ matchCountFunctionList_BASE2.map{ case (id,arr, varFcn,fcn) => {
      ("CPASSN_"+id,Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => {
        varFcn(vc) && vc.filter != "." && vc.filter.split(",")(1) != "." && vc.filter.split(",")(0) != "."
      },fcn)
    }}
    
    val matchCountFunctionList_FINAL : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int,Int) => Boolean)] = matchCountFunctionList_BASE2 ++ matchCountFunctionList_BASE2_CFILTSET
    
    ////////////////////////////////////////////////////////////////////
    
    val mmCountFunctionList_BASE : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)] = Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)](
        ("noCall",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx: Int, i : Int) => {
              ! isCalled(vc,idx,i)
        }),
        ("called",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx: Int, i : Int) => {
              isCalled(vc,idx,i)
        }),
        ("HomRef",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx: Int, i : Int) => {
              isCalled(vc,idx,i) & isRef(vc,idx,i)
        }),
        ("Het",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx: Int, i : Int) => {
              isCalled(vc,idx,i) & isHet(vc,idx,i)
        }),
        ("HomAlt",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx: Int, i : Int) => {
              isCalled(vc,idx,i) & isHomAlt(vc,idx,i)
        }),
        ("Other",Array.fill[Int](matchIdx.length)(0), (vc : SVcfVariantLine) => true,(vc : SVcfVariantLine, idx: Int, i : Int) => {
              isCalled(vc,idx,i) & isOther(vc,idx,i)
        })
    );
    val mmCountFunctionList_BYTYPE  : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)]  = mmCountFunctionList_BASE.map{ case (id,arr,varFcn,fcn) => {
      (id + "_SNV",Array.fill[Int](matchIdx.length)(0),(vc : SVcfVariantLine) => {
        isSNV(vc);
      },fcn)
    }} ++ mmCountFunctionList_BASE.map{ case (id,arr,varFcn,fcn) => {
      (id + "_INDEL",Array.fill[Int](matchIdx.length)(0),(vc : SVcfVariantLine) => {
        ! isSNV(vc);
      },fcn)
    }}
    
    val mmCountFunctionList_BYSWAP  : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)] =  mmCountFunctionList_BASE.flatMap{ case (id,arr,varFcn,fcn) => {
       baseSwapList.map{ case ((r1,a1),(r2,a2)) => {
         (id + "_SNV",Array.fill[Int](matchIdx.length)(0),(vc : SVcfVariantLine) => {
           isSNV(vc) && ((vc.ref == r1 && vc.alt == a1) || (vc.ref == r2 && vc.alt == a2))
         },fcn)
       }}
    }}
    
    
    val mmCountFunctionList_BASE2  : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)] = mmCountFunctionList_BASE ++ mmCountFunctionList_BYTYPE ++ mmCountFunctionList_BYSWAP;
    
    val mmCountFunctionList_CPASS  : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)] = mmCountFunctionList_BASE2.map{ case (id,arr,varFcn,fcn) => {
      ("CPASS_"+id,Array.fill[Int](matchIdx.length)(0),(vc : SVcfVariantLine) => {
        varFcn(vc) && vc.filter == "."
      },fcn)
    }}
    
    val mmCountFunctionList_A  : Vector[(String, Array[Int], SVcfVariantLine => Boolean, (SVcfVariantLine,Int,Int) => Boolean)] = mmCountFunctionList_BASE2 ++ mmCountFunctionList_CPASS
    
    val mmCountFunctionList_B = mmCountFunctionList_A.map{ case (id,arr,varFcn,fcn) => {
      (id,Array.fill[Int](matchIdx.length)(0),varFcn,fcn)
    }}

    ////////////////////////////////////////////////////////////////////

    def iteratePos() : Seq[SVcfVariantLine] = {
      val (chrom,pos) = (vcIter1.head.chrom,vcIter1.head.pos);
      val v1s = extractWhile(vcIter1)( vc => { vc.pos ==  pos && vc.chrom == chrom});
      val v2s = extractWhile(vcIter2)( vc => { vc.pos ==  pos && vc.chrom == chrom});
      val alts : Vector[String] = (v1s.map{vc => vc.alt.head}.toSet ++ v2s.map{vc => vc.alt.head}).toSet.toVector.sorted;
      reportln("      Found "+alts.length+" alt alleles at this position: ["+alts.mkString(",")+"]","debug")
      
      alts.map{ alt => {
        val v1idx = v1s.indexWhere{ vc => { vc.alt.head == alt }};
        val v2idx = v2s.indexWhere{ vc => { vc.alt.head == alt }};
        if(v1idx == -1 ){
          val vc = v2s(v2idx).getOutputLine();
          alleAt2not1ct += 1;
          reportln("      Allele found only on iter2: "+alt,"debug")
          writeVC2(vc);
        } else if(v2idx == -1) {
          val vc = v1s(v1idx).getOutputLine();
          alleAt1not2ct += 1;
          reportln("      Allele found only on iter1: "+alt,"debug")
          writeVC1(vc);
        } else {
          val vc1 = v1s(v1idx).getOutputLine();
          val vc2 = v2s(v2idx).getOutputLine();
          matchCt += 1;
          reportln("      Allele found on both: "+alt,"debug")
          writeJointVC(vc1,vc2);
        }
      }}
    }
    
    def writeVC1(vc : SVcfVariantLine) : SVcfVariantLine = {
      val gt = vc.genotypes;
      val gto = SVcfGenotypeSet(vc.format, gt.genotypeValues.map{ ga => { matchIdx.map{ case (samp,idx1,idx2) => {
        ga(idx1);
      }}}.toArray});
      val vco = SVcfOutputVariantLine(
       in_chrom = vc.chrom,in_pos = vc.pos,in_id = vc.id,in_ref = vc.ref,in_alt = vc.alt,in_qual = vc.qual,in_filter = vc.filter, in_info = vc.info,
       in_format = vc.format,
       in_genotypes = gto
      )
      out1.write(vco.getVcfString+"\n");
      val idx = vco.format.indexOf(genoTag1)
      mmCountFunctionList_A.foreach{ case (id,arr,varFcn,fcn) => {
        if(varFcn(vco)){
          matchIdxIdx.foreach{ i => {
            if(fcn(vc,idx,i)) arr(i) += 1;
          }}
        }
      }}
      vco;
    }
    def writeVC2(vc : SVcfVariantLine) : SVcfVariantLine = {
      val gt = vc.genotypes;
      val fmt = vc.format.map{ f => { file2tag + f }}
      val gto = SVcfGenotypeSet(fmt, gt.genotypeValues.map{ ga => { matchIdx.map{ case (samp,idx1,idx2) => {
        ga(idx2);
      }}}.toArray});
      val vco = SVcfOutputVariantLine(
       in_chrom = vc.chrom,in_pos = vc.pos,in_id = vc.id,in_ref = vc.ref,in_alt = vc.alt,in_qual = vc.qual,in_filter = vc.filter, in_info = vc.info,
       in_format = fmt,
       in_genotypes = gto
      )
      out2.write(vco.getVcfString+"\n");
      val idx = vco.format.indexOf(fGenoTag2);
      mmCountFunctionList_B.foreach{ case (id,arr,varFcn,fcn) => {
        if(varFcn(vco)){
          matchIdxIdx.foreach{ i => {
            if(fcn(vc,idx,i)) arr(i) += 1;
          }}
        }
      }}
      vco;
    }

    def writeJointVC(vc1 : SVcfVariantLine, vc2 : SVcfVariantLine) : SVcfVariantLine = {
      val gt1 = vc1.genotypes;
      val fmt1 = vc1.format//.map{ f => { file2tag + f }}
      val gt2 = vc2.genotypes;
      val fmt2 = vc2.format.map{ f => { file2tag + f }}
      val fmt = fmt1 ++ fmt2;
      val gt = SVcfGenotypeSet(fmt, 
        gt1.genotypeValues.map{ ga => { matchIdx.map{ case (samp,idx1,idx2) => {
          ga(idx1);
        }}}.toArray} ++
        gt2.genotypeValues.map{ ga => { matchIdx.map{ case (samp,idx1,idx2) => {
          ga(idx2);
        }}}.toArray}
      );
      val filt = if(vc1.filter == "." && vc2.filter == ".") { "." } else if(vc1.filter == vc2.filter) {
        vc1.filter + "," + vc2.filter;
      } else {
        vc1.filter + "," + vc2.filter;
      }
      val vco = SVcfOutputVariantLine(
       in_chrom = vc1.chrom,in_pos = vc1.pos,in_id = vc1.id,in_ref = vc1.ref,in_alt = vc1.alt,in_qual = vc1.qual,
       in_filter = filt, in_info = vc1.info,
       in_format = fmt,
       in_genotypes = gt
      )
      outM.write(vco.getVcfString+"\n");
      val (idx1,idx2) = (vco.format.indexOf(genoTag1),vco.format.indexOf(fGenoTag2));
      if(idx1 != -1 && idx2 != -1){
        matchCountFunctionList_FINAL.foreach{ case (id,arr,varFcn,fcn) => {
          if(varFcn(vco)){
            matchIdxIdx.foreach{ i => {
              if(fcn(vco,idx1,idx2,i)) arr(i) += 1;
            }}
          }
        }}
      }

      vco;
    }
    
    def iterateMissingVariant(first : Boolean){
      if(first){
        val vc = vcIter1.next();
        writeVC1(vc);
        at1not2ct += 1;
      } else {
        val vc = vcIter2.next().getOutputLine();
        //vc.genotypes.fmt = vc.genotypes.fmt.map{f => {file2tag + f}}
        //out2.write(vc.getVcfString+"\n");
        writeVC2(vc);
        at2not1ct += 1;
      }
    }
  }
  
  
  
  class redoDBNSFP extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "redoDBNSFP", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new UnaryArgument( name = "singleDbFile",
                                         arg = List("--singleDbFile"), // name of value
                                         argDesc = "NOT CURRENTLY SUPPORTED"+
                                                   "" // description
                                       ) ::
                    new BinaryArgument[String](
                                         name = "dbFileDelim", 
                                         arg = List("--dbFileDelim"), 
                                         valueName = "delim",  
                                         argDesc =  ".",
                                         defaultValue = Some("\t")
                                        ) ::
                    new BinaryArgument[String](
                                         name = "chromStyle", 
                                         arg = List("--chromStyle"), 
                                         valueName = "hg19",  
                                         argDesc =  ".",
                                         defaultValue = Some("hg19")
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryArgument[String](
                                         name = "posFieldTitle", 
                                         arg = List("--posFieldTitle"), 
                                         valueName = "pos",  
                                         argDesc =  ".",
                                         defaultValue = Some("pos(1-coor)")
                                        ) ::
                    new BinaryArgument[String](
                                         name = "chromFieldTitle", 
                                         arg = List("--chromFieldTitle"), 
                                         valueName = "chr",  
                                         argDesc =  ".",
                                         defaultValue = Some("chr")
                                        ) ::     
                    new BinaryArgument[String](
                                         name = "altFieldTitle", 
                                         arg = List("--altFieldTitle"), 
                                         valueName = "alt",  
                                         argDesc =  ".",
                                         defaultValue = Some("alt")
                                        ) ::     
                    new BinaryArgument[String](
                                         name = "tagPrefix", 
                                         arg = List("--tagPrefix"), 
                                         valueName = "prefix",  
                                         argDesc =  ".",
                                         defaultValue = Some("SWH_dbNSFP_")
                                        ) :: 
                    new BinaryOptionArgument[List[String]](
                                         name = "keepTags", 
                                         arg = List("--keepTags"), 
                                         valueName = "tag1,tag2,...",  
                                         argDesc =  ""
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "dropTags", 
                                         arg = List("--dropTags"), 
                                         valueName = "tag1,tag2,...",  
                                         argDesc =  ""
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "dbnsfpfile",
                                         valueName = "dbNSFP3.0b2a.txt.gz",
                                         argDesc = "A gene annotation GTF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         RedoDBNSFPannotation(
             dbnsfpfile = parser.get[String]("dbnsfpfile"),
             chromStyle = parser.get[String]("chromStyle"), 
             chromList = parser.get[Option[List[String]]]("chromList"),
            // filesByChrom = ! parser.get[Boolean]("singleDbFile"),
             posFieldTitle = parser.get[String]("posFieldTitle"),
             chromFieldTitle = parser.get[String]("chromFieldTitle"),
             altFieldTitle = parser.get[String]("altFieldTitle"),
             singleDbFile = parser.get[Boolean]("singleDbFile"),
             dbFileDelim = parser.get[String]("dbFileDelim"),
             tagPrefix = parser.get[String]("tagPrefix"),
             dropTags = parser.get[Option[List[String]]]("dropTags"),
             keepTags = parser.get[Option[List[String]]]("keepTags")
         ).walkVCFFile(
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList")
         )
       }
     }
    
  }

  case class RedoDBNSFPannotation( dbnsfpfile : String, chromStyle : String, chromList : Option[List[String]], 
                                   posFieldTitle : String,chromFieldTitle : String, altFieldTitle : String,
                                   singleDbFile : Boolean,
                                   dbFileDelim : String,
                                   tagPrefix : String,
                                   dropTags : Option[List[String]],keepTags : Option[List[String]]
                                  ) extends internalUtils.VcfTool.VCFWalker {
    
    //if(! filesByChrom){
    //  error("Fatal error: operation with the --singleDbFile flag is not yet supported!");
    //}
    
    var currChrom = "chr20";

    var fileCells = if(singleDbFile){
      getLinesSmartUnzip(dbnsfpfile).map(_.split(dbFileDelim))
    } else {
      getLinesSmartUnzip(dbnsfpfile + currChrom).map(_.split(dbFileDelim))
    }
    val dbHeader = fileCells.next.zipWithIndex.map{case (s,i) => if(s.head == '#') s.tail else s}.map(s => { s.trim() });
    
    reportln("DBNSFP file ("+dbnsfpfile+") header: ","debug");
    dbHeader.zipWithIndex.foreach{ case (s,i) => reportln("    "+i+"=\""+s+"\"","debug")}
    
    val keyMap = dbHeader.zipWithIndex.toMap;
    val altMap = Map[String,Int](("A" -> 0),("C" -> 1),("G" -> 2),("T" -> 3));
    
    val posIdx = keyMap(posFieldTitle);
    val alleIdx = keyMap("alt");
    
    val zeroString = Array.ofDim[String](dbHeader.size).map(_ => ".");
    
    var currPositionMap : Seq[(String,Array[String])] = Seq[(String,Array[String])]()
    
    val chromFieldIdx = keyMap(chromFieldTitle);
    val posFieldIdx = keyMap(posFieldTitle);
    val altFieldIdx = keyMap(altFieldTitle);
    
    var currIterator : BufferedIterator[(String,Int,String,Array[String])] = fileCells.map( cells => {
        val chrom = if(singleDbFile) cells(chromFieldIdx) else currChrom;
        val pos = string2int(cells(posFieldIdx));
        val alle = cells(altFieldIdx);
        (chrom,pos,alle,cells)
      }).buffered
    
    //var currCells = Array.ofDim[String](dbHeader.length) //currReader.next.split("\t");
    //var currPos = string2int(currCells(keyMap("hg19_pos(1-based)")));
    var currPos = -1;
    var lastRequestedPos = -1;
    
    def setPos(){
      if(currIterator.hasNext){
        val (chr, pos, alle, cells) = currIterator.head;
        val currPosVector = extractWhile(currIterator){ case (ch,p,a,c) => { ch == chr && p == pos } }
        
        //val (currPosIter, remainderIter) = currIterator.span{ case (p,a,c) => { p == pos }};
        //currIterator = remainderIter;
        currPositionMap = currPosVector.map{case (chr,p,a,c) => ((a,c))}
        currPos = pos;
        currChrom = chr;
      }
    }
    
    def shiftToPosition(chrom : String, pos : Int) : Boolean = {
      if(chrom == currChrom){
        if(currPos > pos){
          if(lastRequestedPos > pos){
            warning("Illegal backward reference in DBNSFP file! Is VCF sorted? ("+chrom+":"+pos+")","DBNSFP_BACKREF",100);
          }
          lastRequestedPos = pos;
          return false;
        } else if(currPos == pos){
          //do nothing!
          lastRequestedPos = pos;
          return true;
        } else {
          //currIterator = currIterator.dropWhile{ case (p,alle,cells) => p < pos }
          skipWhile(currIterator){ case (chr,p,alle,cells) => chr == chrom && p < pos }
          if(currIterator.head._1 != chrom){
            lastRequestedPos = pos;
            return false;
          }
          setPos();
          lastRequestedPos = pos;
          return pos == currPos && chrom == currChrom;
        }
      } else {
        reportln("Switching to chromosome: "+currChrom +" ["+ getDateAndTimeString+"]","debug");          
        if(singleDbFile){
          skipWhile(currIterator){ case (chr,p,alle,cells) => chr != chrom }
          if(! currIterator.hasNext){
            reportln("Reached file-end, returning to start of file... "+" ["+ getDateAndTimeString+"]","debug");
            currIterator = getLinesSmartUnzip(dbnsfpfile).drop(1).map( line => {
                 val cells = line.split(dbFileDelim);
                 val chrom = cells(chromFieldIdx);
                 val pos = string2int(cells(posFieldIdx));
                 val alle = cells(altFieldIdx);
                 (chrom,pos,alle,cells)
            }).buffered
            skipWhile(currIterator){ case (chr,p,alle,cells) => chr != chrom }
            if(! currIterator.hasNext){
              warning("Chromosome "+chrom+" not found!","CHROM_NOT_FOUND_WARNING",100);
            }
          }
        } else {
          currIterator = getLinesSmartUnzip(dbnsfpfile + chrom).drop(1).map( line => {
            val cells = line.split(dbFileDelim);
            val pos = string2int(cells(keyMap(posFieldTitle)));
            val alle = cells(keyMap("alt"));
            (chrom,pos,alle,cells)
          }).buffered
        }
        currPos = -1;
        currChrom = chrom;
        lastRequestedPos = -1;
        reportln("Switched to chromosome: "+currChrom +" ["+ getDateAndTimeString+"]","debug");          
        return shiftToPosition(chrom,pos);
      }
    }
    
    val tagsToWrite : Seq[String] = dbHeader.filter(t => {
      dropTags match {
        case Some(dt) => {
          ! dt.contains(t);
        }
        case None => true;
      }
    }).filter(t => {
      keepTags match {
        case Some(dt) => {
          dt.contains(t);
        }
        case None => true;
      }
    })
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = { 
      
      val newHeaderLines = List(
        new VCFInfoHeaderLine(tagPrefix+"FoundAtPos", 1, VCFHeaderLineType.Integer, "Equal to 1 if and only if any dbNSFP line was found at the given position."),
        new VCFInfoHeaderLine(tagPrefix+"Found", 1, VCFHeaderLineType.Integer, "Equal to the number of lines found at the given position and matching the given alt allele.")
      ) ++ tagsToWrite.toList.zipWithIndex.map{case (title,i) => {
        new VCFInfoHeaderLine(tagPrefix+title, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Info from column "+title+" (col "+i+") of dbNSFP file ("+dbnsfpfile+")");
      }}
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      
      ((vcIter.map(vc => walkVC(vc)),newHeader));
      
    }
    
    def walkVC(vc : VariantContext) : VariantContext = {
      var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(vc);
      
      val chrom = vc.getContig();
      val pos = vc.getStart();
      val isInDB = shiftToPosition(chrom,pos);
      //val altAlleles = Range(0,vc.getNAlleles()-1).map((a) => vc.getAlternateAllele(a)).zipWithIndex.filter{case (alt,altIdx) => { alt.getBaseString() != "*" }}
      val altAllelesRaw = internalUtils.VcfTool.getAllelesInOrder(vc).toVector.tail.zipWithIndex
      val altAlles = altAllelesRaw.filter{case (alt,altIdx) => { alt.getBaseString() != "*" }}
      if(altAlles.length > 1){
        error("Fatal Error: wrong number of alt alleles!")
      }
      if(isInDB) vb = vb.attribute(tagPrefix+"FoundAtPos","1");
      else vb = vb.attribute(tagPrefix+"FoundAtPos","0");
      
      if(isInDB && altAlles.length > 0 && currPositionMap.exists{ case (key,c) => key == altAlles.head._1.getBaseString() } ){
        val alt = altAlles.head._1.getBaseString();
        val matches = currPositionMap.filter{ case (key,c) => key == alt }.map{ case (key,c) => c};
        
        if(matches.length == 1){
          matches.head.zip(dbHeader).filter{ case (v,title) => {
            tagsToWrite.contains(title)
          }}.foreach{ case (v,title) => {
            vb = vb.attribute(tagPrefix+title,cleanInfoField(v));
          }}
        } else {
          dbHeader.zipWithIndex.filter{ case (tag,idx) => {
            tagsToWrite.contains(tag);
          }}.foreach{ case (tag,idx) => {
            vb = vb.attribute(tagPrefix+tag,matches.map{ case varray => cleanInfoField(varray(idx)) }.filter( v => v != "." ).padTo(1,".").mkString(","));
          }}
        }
 
        vb = vb.attribute(tagPrefix+"Found",matches.length);
      } else {
        tagsToWrite.foreach(title => {
          vb = vb.attribute(tagPrefix+title,".");
        })
        vb = vb.attribute(tagPrefix+"Found","0");
      }
      
      return vb.make();
    }
    
    //def cleanInfoField(s : String) : String = {
    //  var out = s;
    //  
    //  out = out.replaceAll("\\||\\.|;",",").replaceAll("/| |-|:","_").replaceAll("[\\(\\)\\[\\]]","").replaceAll("[_]+","_").replaceAll("[,]+",",").replaceAll("_$|,$","");
    //  
    //  return out;
    //}
    def cleanInfoField(s : String) : String = {
      var out = s;
      
      //out = out.replaceAll("[;|.]+",",").replaceAll("[/ ;-]","_").replaceAll("[\\(\\)\\[\\]]","").replaceAll("[_,]+$|^[_,]+","");
//      out = out.replaceAll("[/ :-]+","_").replaceAll("[\\(\\)\\[\\]]","").split("[;|.]+").map(k => k.replaceAll("[_]+$|^[_]+","")).filter(k => k != "").padTo(1,".").mkString(",");

      out = out.replaceAll("[/ :-]+","_").replaceAll("[\\(\\)\\[\\]]","").split("[;|,]+").map(k => k.replaceAll("[_]+$|^[_]+","")).filter(k => k != "").padTo(1,".").mkString(",");
      return out;
    }

    
  }
  
  
  
  
  class addTXAnno extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "addTxInfoToVCF", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility adds an array of new VCF tags with information about the transcriptional changes caused by each variant. "+BETA_WARNING,
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "inputSavedTxFile", 
                                         arg = List("--inputSavedTxFile"), 
                                         valueName = "txdata.data.txt.gz",  
                                         argDesc =  "Loads a saved TXdata file. Either this parameter OR the --genomeFA parameter must be set. Using this file will be much faster than regenerating the tx data from the gtf/fasta."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "genomeFA", 
                                         arg = List("--genomeFA"), 
                                         valueName = "genome.fa.gz",  
                                         argDesc =  "The genome fasta file. Can be gzipped or in plaintext. Either this parameter OR the --inputSavedTxFile parameter must be set!"
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "summaryFile", 
                                         arg = List("--summaryFile"), 
                                         valueName = "filename.txt",  
                                         argDesc =  "An optional extra output file that contains debugging information."
                                        ) ::
                    new UnaryArgument( name = "cdsRegionContainsStop",
                                         arg = List("--cdsRegionContainsStop"), // name of value
                                         argDesc = "Use this flag if the input GTF annotation file includes the STOP codon in the CDS region. Depending on the source of the annotation file, some GTF files include the STOP codon, some omit it. The UCSC knowngenes annotation file does NOT include CDS regions."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "addSummaryCLNSIG",
                                         arg = List("--addSummaryCLNSIG"), // name of value
                                         argDesc = "Special-purpose flag for use with specialized ClinVar VCFs. NOT FOR GENERAL USE!"+
                                                   "" // description
                                       ) ::
                    new BinaryOptionArgument[String]( name = "addCanonicalTags",
                                         arg = List("--addCanonicalTags"), // name of value
                                         valueName = "knownCanonical.txt",
                                         argDesc = "Supply a list of canonical transcripts, add tags that indicate canonical-transcript-only variant info."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "splitMultiAllelics",
                                         arg = List("--splitMultiAllelics"), // name of value
                                         argDesc = "If this flag is used, multiallelic variants will be split into multiple separate VCF lines."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "geneVariantsOnly",
                                         arg = List("--geneVariantsOnly"), // name of value
                                         argDesc = "If this flag is used, only variants that fall on or near known genes will be written."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "nonNullVariantsOnly",
                                         arg = List("--nonNullVariantsOnly"), // name of value
                                         argDesc = "If this flag is used, only write variants that have non-null alt alleles."+
                                                   "" // description
                                       ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "addBedTags", 
                                         arg = List("--addBedTags"), 
                                         valueName = "TAGTITLE:filedesc:bedfile.bed,TAGTITLE2:filedesc2:bedfile2.bed",  
                                         argDesc =  "List of tags and bed files that define said tags. For each tag, the variant will have a tag value of 1 iff the variant appears on the bed file region, and 0 otherwise."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "txInfoFile", 
                                         arg = List("--txInfoFile"), 
                                         valueName = "txInfoFile.txt",  
                                         argDesc =  "Outputs an optional debugging file."
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

                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "gtffile.gtf.gz",
                                         argDesc = "A gene annotation GTF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
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
                       splitMultiAllelics = parser.get[Boolean]("splitMultiAllelics"),
                       addSummaryCLNSIG = parser.get[Boolean]("addSummaryCLNSIG"),
                       addCanonicalTags = parser.get[Option[String]]("addCanonicalTags"),
                       geneVariantsOnly = parser.get[Boolean]("geneVariantsOnly"),
                       nonNullVariantsOnly = parser.get[Boolean]("nonNullVariantsOnly"),
                       addBedTags = parser.get[Option[List[String]]]("addBedTags")
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
                addSummaryCLNSIG : Boolean,
                addCanonicalTags : Option[String],
                geneVariantsOnly : Boolean,
                nonNullVariantsOnly : Boolean,
                addBedTags : Option[List[String]],
                bufferSize : Int = 32, 
                vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ){
    
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(vcffile, 
                                       chromList = chromList,
                                       vcfCodes = vcfCodes);
    
    val summaryWriter = if(summaryFile.isEmpty) None else Some(openWriterSmart(summaryFile.get));
    
    val walkers = Seq[VCFWalker]( 
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
                geneVariantsOnly = geneVariantsOnly,
                bufferSize = bufferSize,
                addBedTags = addBedTags,
                vcfCodes =vcfCodes
                )
        ) ++ (
            if(nonNullVariantsOnly){
              Seq[VCFWalker]( FilterNonVariantWalker() );
            } else {
              Seq[VCFWalker]();
            }
        ) ++ (
            if(groupFile.isEmpty){
              Seq[VCFWalker]();
            } else {
              Seq[VCFWalker](AddGroupInfoAnno(groupFile = groupFile, groupList = None, superGroupList  = superGroupList, chromList = chromList))
            }
        ) ++ (
            if(splitMultiAllelics){
              Seq[VCFWalker](SplitMultiAllelics(vcfCodes = vcfCodes, clinVarVariants = false, splitSimple = false));
            } else {
              Seq[VCFWalker]();
            }
        ) ++ (
            if(addSummaryCLNSIG){
              Seq[VCFWalker](AddSummaryCln(vcfCodes = vcfCodes));
            } else {
              Seq[VCFWalker]();
            }
        ) ++ (
            if(addCanonicalTags.isDefined){
              Seq[VCFWalker](AddCanonicalInfo(canonicalTxFile = addCanonicalTags.get));
            } else {
              Seq[VCFWalker]();
            }
        )
    
        
    val (finalIter, finalHeader) = walkers.foldLeft((vcIter,vcfHeader)){case ((iter,header),walker) => {
      walker.walkVCF(vcIter = iter, vcfHeader=header)
    }}
        
    val vcfWriter = internalUtils.VcfTool.getVcfWriter(outfile, header = finalHeader);

    finalIter.foreach(vc => {
      vcfWriter.add(vc)
    })
    vcfWriter.close();
    
    if(! summaryWriter.isEmpty) summaryWriter.get.close();
  }
  
  case class FilterNonVariantWalker() extends internalUtils.VcfTool.VCFWalker {
    def walkVCF(vcIter : Iterator[VariantContext],vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      (vcIter.filter{vc => {
        vc.getAlternateAlleles().asScala.filter(_.getBaseString() != "*").length > 0;
      }},vcfHeader)
    }
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
                geneVariantsOnly : Boolean,
                chromList : Option[List[String]],
                txToGeneFile : Option[String],
                bufferSize : Int = 32, 
                addBedTags : Option[List[String]],
                vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ) extends internalUtils.VcfTool.VCFWalker {
      
      reportln("Starting TX parse...","progress");
       
      val bedTags : Seq[(String,String,internalUtils.commonSeqUtils.GenomicInterval => Boolean)] = addBedTags match {
        case Some(bt) => {
          bt.map{b => {
            val pair : Array[String] = b.split(":");
            if(pair.length != 3) error("Each comma-delimited element of parameter addBedTags must have exactly 3 colon-delimited elements (tag:desc:filename.bed).")
            val (t,desc,f) : (String,String,String) = (pair(0),pair(1),pair(2));
            val arr : internalUtils.genomicAnnoUtils.GenomicArrayOfSets[String] = internalUtils.genomicAnnoUtils.GenomicArrayOfSets[String](false);
            val lines = getLinesSmartUnzip(f);
            lines.foreach(line => {
              val cells = line.split("\t");
              val (chrom,start,end) = (cells(0),string2int(cells(1)),string2int(cells(2)))
              arr.addSpan(internalUtils.commonSeqUtils.GenomicInterval(chrom, '.', start,end), "CE");
            })
            arr.finalizeStepVectors;
            val isOnBedFunc : (internalUtils.commonSeqUtils.GenomicInterval => Boolean) = {
              (iv : internalUtils.commonSeqUtils.GenomicInterval) => {
                ! arr.findIntersectingSteps(iv).foldLeft(Set[String]()){case (soFar,(iv,currSet)) => {
                  soFar ++ currSet;
                }}.isEmpty
              }
            }
            (t,desc,isOnBedFunc)
          }}
        }
        case None => {
          Seq();
        }
      }
      
      /*
       * 
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
       * 
       */
      
      
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
      }) ++ bedTags.map{ case (tagString,descString,bedFunction) => {
          new VCFInfoHeaderLine(tagString, 1, VCFHeaderLineType.Integer, "Variant is found on bed file "+descString)
      }}
      
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
      
      
      return ((
        vcIter.map(v => {
             annotateVcfStreamOpt(v,summaryWriter,vcfCodes,bufferSize,txgaos,TXSeq,txToGene,geneVariantsOnly=geneVariantsOnly,bedTags=bedTags)
        }).filter(_.isDefined).map(_.get), newHeader ));
      }
    }
   
  def annotateVcfStreamOpt(v : VariantContext, writer : Option[WriterUtil], vcfCodes : VCFAnnoCodes,
                        bufferSize : Int, txgaos : GenomicArrayOfSets[String],
                        TXSeq : scala.collection.mutable.Map[String,TXUtil],
                        txToGene : Option[(String => String)],
                        geneVariantsOnly : Boolean,
                        bedTags : Seq[(String,String,internalUtils.commonSeqUtils.GenomicInterval => Boolean)]) : Option[VariantContext] = {
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
        if(geneVariantsOnly && txList.length == 0) return None;
        
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
        
        if(! bedTags.isEmpty){
          val variantIV = internalUtils.commonSeqUtils.GenomicInterval(v.getContig(),'.', start = v.getStart() - 1, end = math.max(v.getEnd(),v.getStart+1));
          bedTags.foreach{ case (tagString,desc,bedFunction) => {
            vb = vb.attribute(tagString, if(bedFunction(variantIV)) "1" else "0");
          }}
        }
      }
      return Some(vb.make());
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
    
    writer.write("#chrom\tstart\tend\tprotSymbol\tprotLen\ttxID\tdomainStartAA\tdomainEndAA\tdomainID\ttxLenAA\ttxStart\ttxEnd\tdomainUID\n")
    
    var badSpanSet = Set[(String,Int,Int,String)]();
    var unkTx = Set[String]();
    
    var uidSet = Set[String]();
    
    for((line,lnct) <- lines.zipWithIndex){
      val cells = line.split("\t");
      val txid  = cells(txIdx);
      val protLenAA = string2int(cells(1));
      //val txLenAA = string2int(cells(7));
      val startAA = string2int(cells(startIdx)) - 1;
      val endAA   = if(cells(endIdx).head == '>') string2int(cells(endIdx).tail) else string2int(cells(endIdx))
      val startC = startAA * 3;
      val endC   = endAA * 3;
      
      val cleanID = cells(5).replaceAll("\\||/| |\\.|,|-|:|;","_").replaceAll("[\\(\\)\\[\\]]","").replaceAll("[_]+","_").replaceAll("_$","");
      //cells(5).replaceAllLiterally(" ","_").replaceAllLiterally("|","_").replaceAllLiterally(".","_").replaceAllLiterally("-","_").replaceAllLiterally(",","_").replaceAllLiterally("/","_").replaceAllLiterally("(","").replaceAllLiterally(")","")
      var idNum = 1;
      var uid = cells(0) + ":" + cleanID + ":" + idNum;
      while(uidSet.contains(uid)){
        idNum = idNum + 1;
        uid = cells(0) + ":" + cleanID + ":" + idNum;
      }
      
      uidSet = uidSet + uid;
      
      if(TXSeq.contains(txid)){
          TXSeq.get(txid) match {
            case Some(tx) => {
              val txLenAA = tx.cLen/3;
              if(protLenAA + 1 == txLenAA){
                val gSpans = tx.convertCSpanToGSpan(startC,endC);
                
                gSpans.foreach{case (s,e) => {
                  writer.write(tx.chrom+"\t"+s+"\t"+e+"\t"+line+"\t"+txLenAA+"\t"+tx.gStart+"\t"+tx.gEnd+"\t"+uid+"\n");
                }}
              } else {
                badSpanSet = badSpanSet + ((tx.chrom,tx.gStart,tx.gEnd,txid)); 
              }
            }
            case None => {
              //writer.write(line+"\t-1\t-1\t-1\n");
              //impossible state!
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
        w.write("chrom\tstart\tend\ttxID\n")
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
        (if(addCounts && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpAltAC_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of alleles called for the alt allele(s)")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpAltHetCt_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as heterozygous for the alt allele(s)")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpAltHomCt_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as homozygous for the alt allele(s)")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpRefHomCt_TAG + g, 1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as homozygous for the reference allele")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpRefHetCt_TAG + g, 1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as heterozygous for the reference allele")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpRefAC_TAG + g, 1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of alleles called for the reference allele")) else List[VCFHeaderLine]()) ++
        List[VCFHeaderLine]()
      }).flatten.toList
      
      //        (if(addFreq && addAlle)     List(new VCFInfoHeaderLine(vcfCodes.grpRefAF_TAG + g, 1, VCFHeaderLineType.Float, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the proportion of alleles called for the reference allele")) else List[VCFHeaderLine]()) ++

  
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
              if(sampSet.contains(samp) && g.countAllele(alle) == 1){
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
          
          if(addCounts && addAlle)   vb = vb.attribute(vcfCodes.grpAltAC_TAG + g,    groupAlleCts.tail.map(_(i)).toList.asJava);
          if(addCounts && addHetHom) vb = vb.attribute(vcfCodes.grpAltHetCt_TAG + g, groupHetCts.tail.map(_(i)).toList.asJava);
          if(addCounts && addHetHom) vb = vb.attribute(vcfCodes.grpAltHomCt_TAG + g, groupHomCts.tail.map(_(i)).toList.asJava);
          if(addCounts && addHetHom) vb = vb.attribute(vcfCodes.grpRefHomCt_TAG + g, groupHomCts.head(i).toString());
          if(addCounts && addHetHom) vb = vb.attribute(vcfCodes.grpRefHetCt_TAG + g, groupHetCts.head(i).toString());
          if(addCounts && addAlle) vb = vb.attribute(vcfCodes.grpRefAC_TAG + g, groupAlleCts.head(i).toString());
          
          /*
           (if(addCounts && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpAltAC_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of alleles called for the alt allele(s)")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpAltHetCt_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as heterozygous for the alt allele(s)")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpAltHomCt_TAG + g, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as homozygous for the alt allele(s)")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpRefHomCt_TAG + g, 1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as homozygous for the reference allele")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addHetHom)   List(new VCFInfoHeaderLine(vcfCodes.grpRefHetCt_TAG + g, 1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of genotypes called as heterozygous for the reference allele")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpRefAC_TAG + g, 1, VCFHeaderLineType.Integer, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the number of alleles called for the reference allele")) else List[VCFHeaderLine]()) ++
        (if(addCounts && addAlle)   List(new VCFInfoHeaderLine(vcfCodes.grpRefAF_TAG + g, 1, VCFHeaderLineType.Float, "For the "+groupToSampleMap(g).size+" samples in group "+g+", the proportion of alleles called for the reference allele")) else List[VCFHeaderLine]()) ++
        
           */
          
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
  
  
  val NUMREPORTBADLEN=100;
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
              val atr = attr.split(",",-1);
              if(atr.length != numAlle - 1){
                warning("WARNING: ATTR: \""+key+"\"=\""+attr+"\" of type \"A\"\n   VAR:\n      atr.length() = "+atr.length + " numAlle = "+numAlle+
                    (if(internalUtils.optionHolder.OPTION_DEBUGMODE) "\n   "+vc.toStringWithoutGenotypes() else ""),
                    "INFO_LENGTH_WRONG",NUMREPORTBADLEN);
                repToSeq(atr.mkString(","),numAlle - 1);
              } else {
                atr.toSeq;
              }
            }
            (key,a)
          }}
          val attrR = copyAttrR.map{case (key,attr) => {
            val a = {
              val atr = attr.split(",",-1);
              if(atr.length != numAlle){
                warning("WARNING: ATTR: \""+key+"\"=\""+attr+"\" of type \"R\"\n   VAR:\n      atr.length() = "+atr.length + " numAlle = "+numAlle+
                       (if(internalUtils.optionHolder.OPTION_DEBUGMODE) "\n   "+vc.toStringWithoutGenotypes() else ""),
                       "INFO_LENGTH_WRONG",NUMREPORTBADLEN);
                repToSeq(atr.mkString(","),numAlle);
              } else {
                atr.toSeq;
              }
            }
            (key,a)
          }}

          val ANNSTR =  trimBrackets(attr.getOrElse("ANN","")).split(",",-1).map(ann => { (ann.trim.split("\\|",-1)(0),ann.trim) });
          
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
            
            var malformedADCT = 0;
            
            if(vc.hasGenotypes()){
            vb.genotypes(GenotypesContext.create(
                new java.util.ArrayList[Genotype](vc.getGenotypesOrderedBy(sampNames).asScala.map{gt => {
                  var gb = new GenotypeBuilder(gt.getSampleName());
                  if(gt.hasAD()){
                    val ad = gt.getAD();
                    if(ad.length != altAlleles.length + 1){
                      if( malformedADCT < 10 ){
                        warning("Malformed AD tag! AD is: [\""+ad.mkString("\",\"")+"\"], Alles are: [\""+refAlle.getBaseString()+"\"][\""+altAlleles.map{a => a.getBaseString()}.mkString("\",\"")+"\"]\n"+
                                "Original VCF String: "+vc.toStringWithoutGenotypes() +"","Malformed_AD_Tag",200);
                        malformedADCT = malformedADCT + 1;
                      }
                    } else {
                      val refAD = ad.head;
                      val altAD = ad(altAlleIdx);
                      val othAD = ad.zipWithIndex.filter{case (a,i) => {i != 0 && i != altAlleIdx}}.map(_._1).sum
                      gb = gb.AD(Array(refAD,altAD,othAD));
                    }
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
                  val attrCells = attr.split(",",-1);
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
      vb = vb.attribute("CLNERROR","OK");
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
        vb = vb.attribute("ORIG"+key, attr.mkString(","));
        if(attr.length != clnAlle.length){
          isErr = true;
           vb = vb.attribute("CLNERROR","ERR");
        } else {
          val vcfAllelesNotFound = clnMap.count(_ == -1);
          val clnAllelesNotFound = clnAlle.count(i => { i == -1 });
          
          if(vcfAllelesNotFound > 0 && clnAllelesNotFound > 0){
              isErr = true;
              vb = vb.attribute("CLNERROR","ERR");
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
          vb = vb.attribute(key, newAttr.mkString(","));
        }
      })
      vb = vb.attribute("ORIGCLNALLE",v.getAttributeAsList("CLNALLE").asScala.map(_.toString()).mkString(","));
      vb = vb.attribute("CLNALLE",Range(0,altAlleles.length).map(x => (x + 1).toString()).mkString(","));
      if(isErr){
        badCtCt(altAlleles.length) = badCtCt(altAlleles.length) + 1;
      }
      
      return vb.make();
    }
    
    reportln("recodeClinVarCLN() Created...","note")
    
  }
  
  
  case class AddSummaryCln(
                 vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ) extends internalUtils.VcfTool.VCFWalker {
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      
      val newHeaderLines = List(
            new VCFInfoHeaderLine(vcfCodes.CLNVAR_SUMSIG,     1, VCFHeaderLineType.Integer, "Summary clinical significance level, from CLNSIG tag. Collapsed multiple reports into a single reported significance level.")      
      );
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      reportln("Walking input VCF...","note")

      return ((vcIter.map(v => {
         var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
         //val clnsig = v.getAttributeAsList("CLNSIG").asScala.map(_.toString()).toSeq;
         val clnsig = getAttributeAsStringList(v,"CLNSIG")
         
         val cs = internalUtils.CalcACMGVar.getSummaryClinSig(clnsig);
         vb = vb.attribute(vcfCodes.CLNVAR_SUMSIG,cs);
         
         vb.make();
      }), newHeader));
      
    }
  }
  
  

  class AddVariantDomainUtil extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "AddVariantDomains", 
          quickSynopsis = "", 
          synopsis = "", 
          description = " " + ALPHA_WARNING,   
          argList = 
                    new UnaryArgument( name = "countClinVar",
                                         arg = List("--countClinVar"), // name of value
                                         argDesc = "If this flag is used..."+
                                                   "" // description
                                       ) ::       
                    new BinaryOptionArgument[String](
                                         name = "txToGeneFile", 
                                         arg = List("--txToGeneFile"), 
                                         valueName = "txToGene.txt",  
                                         argDesc =  ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "gtfFile", 
                                         arg = List("--gtfFile"), 
                                         valueName = "anno.gtf.gz",  
                                         argDesc =  ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chromList.txt",  
                                         argDesc =  ""
                                        ) ::   
                    new FinalArgument[String](
                                         name = "invcf",
                                         valueName = "variants.vcf",
                                         argDesc = "infput VCF file" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "domainFile",
                                         valueName = "domainFile.txt",
                                         argDesc = "infput domain info file." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "domainSummaryOutFile",
                                         valueName = "domainSummaryOutFile.txt",
                                         argDesc = "output domain info file." // description
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
       
       VcfAnnotateTX.AddVariantDomains(
                   domainFile = parser.get[String]("domainFile"),
                   countClinVar = parser.get[Boolean]("countClinVar"),
                   summaryFile = parser.get[String]("domainSummaryOutFile"),
                   gtfFile = parser.get[Option[String]]("gtfFile"),
                   txToGeneFile = parser.get[Option[String]]("txToGeneFile")
                   ).walkVCFFile(
                   infile    = parser.get[String]("invcf"),
                   outfile   = parser.get[String]("outvcf"),
                   chromList = parser.get[Option[List[String]]]("chromList")
                   )
     }   
    }
  }
  
  case class AddVariantDomains(
                 domainFile : String,
                 gtfFile : Option[String],
                 txToGeneFile : Option[String],
                 countClinVar : Boolean,
                 summaryFile : String,
                 vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                ) extends internalUtils.VcfTool.VCFWalker {
    
    
    reportln("Reading txToGene map... ["+internalUtils.stdUtils.getDateAndTimeString+"]","debug");
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
    reportln("Reading geneToTX map... ["+internalUtils.stdUtils.getDateAndTimeString+"]","debug");
    
    reportln("Reading txdata ... ["+internalUtils.stdUtils.getDateAndTimeString+"]","debug");
    val txData = gtfFile match {
      case Some(txf) => {
          val ipr = internalUtils.stdUtils.IteratorProgressReporter_ThreeLevel("lines", 200, 1000 , 2000 )
          val wrappedIter = internalUtils.stdUtils.wrapIteratorWithProgressReporter(getLinesSmartUnzip(txf) , ipr )
          val txArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](false);
          wrappedIter.foreach{line => {
            val tx = buildTXUtilFromString(line);
            tx.gSpans.foreach{ case (start,end) => {
              txArray.addSpan(GenomicInterval(chromName = tx.chrom, strand = '.', start = start, end = end), tx.txID);
            }}
          }}
         txArray.finalizeStepVectors;
      }
      case None => GenomicArrayOfSets[String](false).finalizeStepVectors;
    }
    reportln("Reading txdata... ["+internalUtils.stdUtils.getDateAndTimeString+"]","debug"); 
    
    val domainClinCount = scala.collection.mutable.AnyRefMap[String,Array[Int]]();
    val domainClinCountBySev = scala.collection.mutable.AnyRefMap[String,Array[Array[Int]]]();
    
    val domainInfo = scala.collection.mutable.AnyRefMap[String,Seq[String]]();
    val domainLoci = scala.collection.mutable.AnyRefMap[String,(String,String,Set[(Int,Int)],(Int,Int))]();
    
    
    val domainTx = scala.collection.mutable.AnyRefMap[String,Set[String]]().withDefault(s => Set[String]());
    val geneArrayTemp : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](false);

    getTableFromLines(getLinesSmartUnzip(domainFile),colNames = Seq("chrom","start","end","protSymbol","domainUID","txID","domainStartAA","domainEndAA")).zipWithIndex.foreach{
      //case (chrom :: start :: end :: protSymbol :: domainUID :: emptyList, lnct) => {
      case (lst,lnct) => {
        val chrom = lst(0);
        val start = string2int(lst(1));
        val end = string2int(lst(2));
        val protSymbol = lst(3);
        val domainUID = lst(4);
        val tx = lst(5);
        val aaStart = lst(6);
        val aaEndRaw = lst(7);
        
        val aaEnd = if(aaEndRaw.head == '>'){
          aaEndRaw.tail;
        } else {
          aaEndRaw;
        }
        
        val iv = GenomicInterval(chromName = chrom, strand = '.', start = start, end = end);
        geneArrayTemp.addSpan(iv,domainUID);
        domainClinCount(domainUID) = Array.ofDim[Int](9);
        domainClinCountBySev(domainUID) = Array.ofDim[Int](2,9);
        domainInfo(domainUID) = Vector(chrom , start.toString() , end.toString());
        if(domainLoci.contains(domainUID)){
          domainLoci(domainUID) = (chrom,tx, domainLoci(domainUID)._3 + ((start,end)), (string2int(aaStart),string2int(aaEnd))   );
        } else {
          domainLoci(domainUID) = (chrom,tx, Set[(Int,Int)]((start,end)), (string2int(aaStart),string2int(aaEnd))   );
        }
        
        val domainTxSet = txData.findIntersectingSteps(iv).foldLeft(Set[String]()){case (soFar,(iv,stepSet)) => soFar ++ stepSet};
        domainTx(domainUID) = domainTx(domainUID) ++ domainTxSet;
      }
      //case (x,lnct) => {
      //  error("Malformed table: "+domainFile + " line "+lnct+" has wrong number of columns ("+x.length+", should be 5)");
      //}
    }
    val geneArray = geneArrayTemp.finalizeStepVectors;
    
    val lvlList = Seq("SYNON","PSYNON","UNK","NONSYNON","PLOF","LLOF");
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      
      val newHeaderLines = List(
            new VCFInfoHeaderLine(vcfCodes.domainIds,      VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "List of domains that the variant intersects with.")      
      );
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      reportln("Walking input VCF...","note")

      return ((addIteratorCloseAction(vcIter.map(v => {
         var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
         val start = v.getStart();
         val chrom = v.getContig();
         val domainList = geneArray.findIntersectingSteps(GenomicInterval(chromName = chrom, strand = '.', start = start, end = start + 1)).foldLeft(Set[String]()){ case (soFar, (iv,currSet)) => {
           soFar ++ currSet;
         }}.toVector.sorted;
         
         if(countClinVar){
           val clnSigRaw = v.getAttributeAsInt(vcfCodes.CLNVAR_SUMSIG, 0);
           val clnSig = if(clnSigRaw == 255) 8 else clnSigRaw;
           domainList.foreach(domainUID => {
             domainClinCount(domainUID)(clnSig) = domainClinCount(domainUID)(clnSig) + 1;
           })
           val varLevel = v.getAttributeAsString(vcfCodes.vMutLVL_TAG,"?");
           val nonSynonStatus = if(varLevel == "NONSYNON" || varLevel == "PLOF" || varLevel == "LLOF") 1 else 0;
           domainList.foreach(domainUID => {
             domainClinCountBySev(domainUID)(nonSynonStatus)(clnSig) = domainClinCountBySev(domainUID)(nonSynonStatus)(clnSig) + 1;
           })
         }
         
         vb = vb.attribute( vcfCodes.domainIds, domainList.mkString(",") );
         //val cs = internalUtils.CalcACMGVar.getSummaryClinSig(clnsig);
         //vb = vb.attribute(vcfCodes.CLNVAR_SUMSIG,cs);
         
         vb.make();
      }), closeAction = () => {
        val writer = openWriterSmart(summaryFile + ".txt.gz");
        val writer2 = openWriterSmart(summaryFile + ".onGene.txt.gz");
        //writer.write("");
        
        val titleLine = "domainUID\tchrom\tgStart\tgEnd\tgenomicSpans\ttxID\taaStart\taaEnd\t"+
                     "txList\tgeneList\t"+
                     Range(0,8).map("ClnSig_"+_).mkString("\t")+"\t"+"ClnSig_255"+"\t"+
                     Range(0,8).map("SYN_ClnSig_"+_).mkString("\t")+"\t"+"SYN_ClnSig_255"+"\t"+
                     Range(0,8).map("NS_ClnSig_"+_).mkString("\t")+"\t"+"NS_ClnSig_255"+"\n";
        writer.write(titleLine);
        writer2.write(titleLine);
        
        domainClinCount.keys.toVector.sorted.foreach(d => {
          val txList = domainTx(d).toVector.sorted;
          val geneList = domainTx(d).map(tx => txToGene(tx)).toSet.toVector.sorted;
            writer.write(d + "\t"+ 
                       domainLoci(d)._1 + "\t" +
                       domainLoci(d)._3.map(_._1).min + "\t"+
                       domainLoci(d)._3.map(_._2).max + "\t"+
                       domainLoci(d)._3.toVector.sorted.map{case (s,e) => { s + "-" + e}}.mkString(",")+"\t"+
                       domainLoci(d)._2 + "\t"+
                       domainLoci(d)._4._1+"\t"+
                       domainLoci(d)._4._2+"\t"+
                       txList.mkString(",") + "\t" + 
                       geneList.mkString(",") + "\t"+
                       domainClinCount(d).mkString("\t")+"\t"+
                       domainClinCountBySev(d)(0).mkString("\t")+"\t"+
                       domainClinCountBySev(d)(1).mkString("\t")+
                       "\n");
        })
        writer.close();

        domainClinCount.keys.toVector.map(d => {
          val geneList = domainTx(d).map(tx => txToGene(tx)).toSet.toVector.sorted;
          val (start,end) = (domainLoci(d)._3.map(_._1).min,domainLoci(d)._3.map(_._2).max);
          (geneList,(start,end),d);
        }).filter{case (geneList,(start,end),d) => {
          geneList.size > 0;
        }}.sortBy{case (geneList,(start,end),d) => {
          (geneList.mkString(","),start,end);
        }}.foreach{ case (geneList,(start,end),d) => {
            val txList = domainTx(d).toVector.sorted;
            writer2.write(d + "\t"+
                       domainLoci(d)._1 + "\t" +
                       domainLoci(d)._3.map(_._1).min + "\t"+
                       domainLoci(d)._3.map(_._2).max + "\t"+
                       domainLoci(d)._3.toVector.sorted.map{case (s,e) => { s + "-" + e}}.mkString(",")+"\t"+
                       domainLoci(d)._2 + "\t"+
                       domainLoci(d)._4._1+"\t"+
                       domainLoci(d)._4._2+"\t"+
                       txList.mkString(",") + "\t" + 
                       geneList.mkString(",") + "\t"+
                       domainClinCount(d).mkString("\t")+"\t"+
                       domainClinCountBySev(d)(0).mkString("\t")+"\t"+
                       domainClinCountBySev(d)(1).mkString("\t")+
                       "\n");
        }}
        writer2.close();
      }), newHeader));
    }
  }
  
  val DEFAULT_CMDADDCANONICALINFO_TAGLIST = DefaultVCFAnnoCodes.txList_TAG + ","+
                                               DefaultVCFAnnoCodes.vMutC_TAG + ","+
                                               DefaultVCFAnnoCodes.vMutP_TAG + ","+
                                               DefaultVCFAnnoCodes.vMutINFO_TAG + ","+
                                               DefaultVCFAnnoCodes.vType_TAG + ","+
                                               DefaultVCFAnnoCodes.vMutR_TAG + ","+
                                               DefaultVCFAnnoCodes.geneIDs;
  val DEFAULT_CMDADDCANONICALINFO_TXTAG = DefaultVCFAnnoCodes.txList_TAG
  
  
  class CmdAddCanonicalInfo extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "AddCanonicalInfo", 
          quickSynopsis = "", 
          synopsis = "", 
          description = " " + ALPHA_WARNING,   
          argList = 
            
            /*
                     vMutC_TAG : String = TOP_LEVEL_VCF_TAG+"varC",
        vMutP_TAG : String = TOP_LEVEL_VCF_TAG+"varPredP",
             */
            
                    new BinaryArgument[String](name = "tagList",
                                           arg = List("--tagList"),  
                                           valueName = "tag1,tag2,tag3,...", 
                                           argDesc = "", 
                                           defaultValue = Some(
                                               DEFAULT_CMDADDCANONICALINFO_TAGLIST
                                           )) ::
                    new BinaryOptionArgument[String](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chromList.txt",  
                                         argDesc =  ""
                                        ) ::   
                    new BinaryArgument[String](name = "txListTag",
                                           arg = List("--txListTag"),  
                                           valueName = "SWH_TXLIST", 
                                           argDesc = "", 
                                           defaultValue = Some(
                                               DEFAULT_CMDADDCANONICALINFO_TXTAG
                                           )) ::
                    new FinalArgument[String](
                                         name = "invcf",
                                         valueName = "variants.vcf",
                                         argDesc = "infput VCF file" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "canonicalTxFile",
                                         valueName = "canonicalTxFile.txt",
                                         argDesc = "Canonical tx file." // description
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
       
       VcfAnnotateTX.AddCanonicalInfo(
                   canonicalTxFile = parser.get[String]("canonicalTxFile"),
                   tagList = parser.get[String]("tagList"),
                   txTag = parser.get[String]("txListTag")
                   ).walkVCFFile(
                   infile    = parser.get[String]("invcf"),
                   outfile   = parser.get[String]("outvcf"),
                   chromList = parser.get[Option[List[String]]]("chromList")
                   )
     }   
    }
  }
  
  case class AddCanonicalInfo(canonicalTxFile : String, tagList : String = DEFAULT_CMDADDCANONICALINFO_TAGLIST, txTag : String = DEFAULT_CMDADDCANONICALINFO_TXTAG) extends internalUtils.VcfTool.VCFWalker {
    val tagSet = tagList.split(",").toSet;
    
    
    val isCanon : (String => Boolean) = {
        val lines = getLinesSmartUnzip(canonicalTxFile);
        val table = getTableFromLines(lines,colNames = Seq("transcript"), errorName = "File "+canonicalTxFile);
        var refSeqSet : Set[String] = table.map(tableCells => {
          val tx = tableCells(0);
          tx
        }).toSet
        reportln("   found: "+refSeqSet.size+" RefSeq transcripts.","debug");
        ((s : String) => {
          refSeqSet.contains(s);
        })
      }
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      val oldHeaderLines = vcfHeader.getInfoHeaderLines().asScala.filter(hl => tagSet.contains(hl.getID()));
      
      val addedHeaderLines = oldHeaderLines.map{ h => {
        if(h.getCountType() == VCFHeaderLineCount.A){
          new VCFInfoHeaderLine(h.getID() + "_CANON", VCFHeaderLineCount.A, h.getType(), "(For the canonical transcript(s) only) "+h.getDescription());
        } else if(h.getCountType() == VCFHeaderLineCount.UNBOUNDED){
          new VCFInfoHeaderLine(h.getID() + "_CANON", VCFHeaderLineCount.UNBOUNDED, h.getType(), "(For the canonical transcript(s) only) "+h.getDescription());
        } else {
          error("FATAL ERROR: ");
          null;
        }
      }}
      val isByAllele = oldHeaderLines.map{ h => {
        if(h.getCountType() == VCFHeaderLineCount.A){
          (h.getID(),true)
        } else if(h.getCountType() == VCFHeaderLineCount.UNBOUNDED){
          (h.getID(),false)
        } else {
          error("FATAL ERROR: ");
          null;
        }
      }}.filter(_._1 != txTag).toVector
      
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader, addedHeaderLines.toVector);
      
      return (vcIter.map{ vc => {
        walkLine(vc,isByAllele=isByAllele,txTag=txTag);
      }}, newHeader);
    }
    
    def walkLine(vc : VariantContext, isByAllele : Vector[(String,Boolean)], txTag : String) : VariantContext = {
      var vb = new VariantContextBuilder(vc);
      
      val txList = getAttributeAsStringList(vc,txTag);
      vb = vb.attribute(txTag + "_CANON",txList.filter(isCanon(_)).padTo(1,".").mkString(","));
      
      if(txList.exists(tx => isCanon(tx))){
        isByAllele.foreach{ case (tag,isA) => {
          val attr = getAttributeAsStringList(vc,tag);
          if(isA){
            vb = vb.attribute(tag+"_CANON",attr.map{a => {
              a.split("\\|").zip(txList).filter{ case (aa,tx) => {
                isCanon(tx);
              }}.map(_._1).mkString("|");
            }}.mkString(","))
          } else {
            vb = vb.attribute(tag+"_CANON",attr.zip(txList).filter{ case (aa,tx) => {
                     isCanon(tx)
                   }}.map(_._1).mkString(","))
          }
        }}
        
      }
      
      return vb.make();
    }
  }
  

  class CmdSummarizeGenotypeStats extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "summarizeGenotypeStats", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new BinaryOptionArgument[Int](
                                         name = "numLinesRead", 
                                         arg = List("--numLinesRead"), 
                                         valueName = "N",  
                                         argDesc =  "Limit file read to the first N lines"
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "master input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[List[String]](
                                         name = "summarizeStatList",
                                         valueName = "stat1,stat2,...",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         SummarizeGenotypeStats(
             summarizeStatList = parser.get[List[String]]("summarizeStatList")
         ).walkVCFFile( 
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList"),
             numLinesRead = parser.get[Option[Int]]("numLinesRead")
         ) 
       }
     }
  }
  
  case class SummarizeGenotypeStats(summarizeStatList : List[String]) extends SVcfWalker {
    val statCommandPairs = summarizeStatList.map(a => {
      val arr = a.split(":"); 
      (arr(0),arr(1))
    })
    
    
    def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader) = {
      //vcfHeader.formatLines = vcfHeader.formatLines ++ Some(new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = filterTag, Number = "1", Type = "Integer", desc = "Equal to 1 if and only if the genotype was filtered due to post-caller quality filters. ("+filter+")"))
      //vcfHeader.addFormatLine(new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = filterTag, Number = "1", Type = "Integer", desc = "Equal to 1 if and only if the genotype was filtered due to post-caller quality filters. ("+filter+")") );
      //vcfHeader.addFormatLine(new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = rawGtTag, Number = "1", Type = "String", desc = "The original GT genotype, prior to genotype-level filtering.") );

      val outIter = vcIter.map{ vc => {
        vc
      }}
      
      return (outIter,vcfHeader);
    }
  }
  
  
  
  
  class CmdFilterGenotypesByStat extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "filterGenotypesByStat", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new BinaryOptionArgument[Int](
                                         name = "numLinesRead", 
                                         arg = List("--numLinesRead"), 
                                         valueName = "N",  
                                         argDesc =  "Limit file read to the first N lines"
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "master input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "filter",
                                         valueName = "filterExpr",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         FilterGenotypesByStat(
             filter = parser.get[String]("filter")
         ).walkVCFFile( 
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList"),
             numLinesRead = parser.get[Option[Int]]("numLinesRead")
         )
       }
     }
  }
  case class FilterGenotypesByStat(filter : String, filterTag : String = "SWHGTFILT", rawGtTag : String = "RAW_GT") extends SVcfWalker {
    val logicParser : internalUtils.VcfTool.SFilterLogicParser[(SVcfVariantLine,Int)] = internalUtils.VcfTool.SGenotypeFilterLogicParser();
    val filterLogic = logicParser.parseString(filter);
    
    def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader) = {
      //vcfHeader.formatLines = vcfHeader.formatLines ++ Some(new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = filterTag, Number = "1", Type = "Integer", desc = "Equal to 1 if and only if the genotype was filtered due to post-caller quality filters. ("+filter+")"))
      vcfHeader.addFormatLine(new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = filterTag, Number = "1", Type = "Integer", desc = "Equal to 1 if and only if the genotype was filtered due to post-caller quality filters. ("+filter+")") );
      vcfHeader.addFormatLine(new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = rawGtTag, Number = "1", Type = "String", desc = "The original GT genotype, prior to genotype-level filtering.") );

      val outIter = vcIter.map{ vc => {
        val vb = vc.getOutputLine();
        //vb.genotypes.genotypeValues = vb.genotypes.genotypeValues :+ Array.fill[String](vb.genotypes.genotypeValues(0).length)("0");
        val ploidy = vb.genotypes.getPloidy();
        val missingGeno = Range(0,ploidy).map{_ => "."}.mkString("/");
        if(! vb.genotypes.fmt.contains(rawGtTag)){
          vb.genotypes.addGenotypeArray( rawGtTag,vb.genotypes.genotypeValues(0).clone() );
        }
        var ftIdx = vb.genotypes.fmt.indexOf(filterTag);
        if(ftIdx == -1){
          vb.genotypes.addGenotypeArray(filterTag,Array.fill[String](vb.genotypes.genotypeValues(0).length)("0"));
          ftIdx = vb.genotypes.fmt.length - 1;
        } //else do nothing. Unfiltered GT's will be overwritten.
        
        vb.genotypes.genotypeValues(0).indices.foreach{ (i) => {
          if(! filterLogic.keep((vc,i))){
            vb.genotypes.genotypeValues(0)(i) = missingGeno;
            vb.genotypes.genotypeValues(ftIdx)(i) = "1";
          }
        }}
        vb;
      }}
      
      return (outIter,vcfHeader);
    }
  }

  case class CalcVariantCountSummary(genotypeTag : String = "GT", keepAltChrom : Boolean = false, singletonGTfile : Option[String] = None,
                                     subFilterExpressionSets : Option[String] = None) {
    
    case class VariantCountSet(
        ctHet : Array[Long],
        ctRef : Array[Long],
        ctAlt : Array[Long],
        ctOthNoAlt : Array[Long],
        ctOthWithAlt : Array[Long],
        ctMis : Array[Long]){
      def addGT(sampGT : String, sampIdx : Int){
          if(sampGT == "./." || sampGT == "."){
            ctMis(sampIdx) = ctMis(sampIdx) + 1;
          } else if(sampGT == "0/0"){
            ctRef(sampIdx) = ctRef(sampIdx) + 1;
          } else if(sampGT == "0/1"){
            ctHet(sampIdx) = ctHet(sampIdx) + 1;
          } else if(sampGT == "1/1"){
            ctAlt(sampIdx) = ctAlt(sampIdx) + 1;
          } else if(sampGT.contains('1')){
            ctOthWithAlt(sampIdx) = ctOthWithAlt(sampIdx) + 1;
          } else {
            ctOthNoAlt(sampIdx) = ctOthNoAlt(sampIdx) + 1;
          }
      }
      def addVC(gt : Array[String]){
        Range(0,gt.length).foreach(idx => addGT(gt(idx),idx));
      }
      def getAll(idx : Int) : String = {
        ctRef(idx) + "\t"+ ctHet(idx) + "\t"+ctAlt(idx)+"\t"+ctOthWithAlt(idx)+"\t"+ctOthNoAlt(idx)
      }
      def getAlt(idx : Int) : String = {
        (ctHet(idx)+ctAlt(idx)+ctOthWithAlt(idx))+"\t"+ctHet(idx) + "\t"+ctAlt(idx)+"\t"+ctOthWithAlt(idx)
      }
      def getNcalled(idx : Int) : Long = {
        ctRef(idx) + ctHet(idx)+ctAlt(idx)+ctOthWithAlt(idx)+ctOthNoAlt(idx)
      }
    }
    
    def makeVariantCountSet(ct : Int) : VariantCountSet = {
      VariantCountSet(
          Array.fill[Long](ct)(0),
          Array.fill[Long](ct)(0),
          Array.fill[Long](ct)(0),
          Array.fill[Long](ct)(0),
          Array.fill[Long](ct)(0),
          Array.fill[Long](ct)(0)
      )
    }
    
    case class VariantCountSetSet(varCt : VariantCountSet,
                                  singCt : VariantCountSet,
                                  indelCt : VariantCountSet,
                                  indelSingCt : VariantCountSet,
                                  snvCt : VariantCountSet,
                                  snvSingCt : VariantCountSet
                                 ){
      def addVariant(gt : Array[String],isSingleton : Boolean, isIndel : Boolean){
        varCt.addVC(gt);
        if(isSingleton){
          singCt.addVC(gt);
        }
        if(isIndel){
          indelCt.addVC(gt);
          if(isSingleton){
            indelSingCt.addVC(gt);
          }
        } else {
          snvCt.addVC(gt);
          if(isSingleton){
            snvSingCt.addVC(gt);
          }
        }
      }
    }
    def makeVariantCountSetSet(ct : Int) : VariantCountSetSet = {
      VariantCountSetSet(varCt = makeVariantCountSet(ct),
                                  singCt  = makeVariantCountSet(ct),
                                  indelCt  = makeVariantCountSet(ct), 
                                  indelSingCt  = makeVariantCountSet(ct),
                                  snvCt  = makeVariantCountSet(ct),
                                  snvSingCt  = makeVariantCountSet(ct)
                                 )
    }
    
    val singletonWriter = singletonGTfile match {
      case Some(f) => {
        Some(openWriterSmart(f))
      }
      case None => None;
    }
    
    def writeSingGT(s : String){
      singletonWriter match {
        case None => {
          //do nothing
        }
        case Some(w) => {
          w.write(s);
        }
      }
    }
    
    def walkVCFFile(infile :String, inputFileList : Boolean, outfile : String, chromList : Option[List[String]], numLinesRead : Option[Int]){
      val (vcIter, vcfHeader) = getSVcfIterators(infile,chromList,numLinesRead,inputFileList);
      val sampleCt = vcfHeader.sampleCt;
      
      val varCt = makeVariantCountSet(sampleCt);
      val singCt = makeVariantCountSet(sampleCt);
      val indelCt = makeVariantCountSet(sampleCt);
      val indelSingCt = makeVariantCountSet(sampleCt);
      val snvCt = makeVariantCountSet(sampleCt);
      val snvSingCt = makeVariantCountSet(sampleCt);
      
      var numVariants = 0;
      var numSingletons = 0;
      var numIndel = 0;
      var numSnv = 0;
      var numFiltVar = 0;
      
      
      
      val subsetList : Seq[(String,SFilterLogic[SVcfVariantLine],VariantCountSetSet)] = subFilterExpressionSets match { 
        case Some(sfes) => {
          val sfesSeq = sfes.split(",")
          sfesSeq.map{ sfe =>
            val sfecells = sfe.split("=").map(_.trim());
            if(sfecells.length != 2) error("ERROR: subfilterExpression must have format: subfiltertitle=subfilterexpression");
            val (sfName,sfString) = (sfecells(0),sfecells(1));
            val parser : SVcfFilterLogicParser = internalUtils.VcfTool.SVcfFilterLogicParser();
            val filter : SFilterLogic[SVcfVariantLine] = parser.parseString(sfString);
            (sfName,filter,makeVariantCountSetSet(sampleCt))
          }
        }
        case None => {
          Seq();
        }
      }
        
        /*subfilterExpression match {
        case Some(filterExpr) => {
          val parser : SVcfFilterLogicParser = internalUtils.VcfTool.SVcfFilterLogicParser();
          val filter : SFilterLogic[SVcfVariantLine] = parser.parseString(filterExpr);
          Some((filter,makeVariantCountSetSet(sampleCt)))
        }
        case None => None;
      }*/
      
      vcIter.filter{vc => vc.genotypes.fmt.contains(genotypeTag) & vc.alt.length > 0}.foreach{vc => {
        if(vc.alt.length > 2) error("Fatal error: multiallelic variant! This utility is only intended to work after splitting multiallelic variants!")
        val alt : String = vc.alt.head;
        val ref : String = vc.ref;
        val gt : Array[String] = vc.getGt(genotypeTag);
        val numAlt = gt.count{ gt => gt.contains('1') }
        val isSingleton = (numAlt == 1);
        val isIndel = ref.length != 1 || alt.length != 1;
        numVariants = numVariants + 1;
          varCt.addVC(gt);
          
          if(isSingleton){
            numSingletons = numSingletons + 1;
            singCt.addVC(gt);
          }
          if(isIndel){
            numIndel = numIndel + 1;
            indelCt.addVC(gt);
            if(isSingleton){
              indelSingCt.addVC(gt);
            }
          } else {
            numSnv = numSnv + 1;
            snvCt.addVC(gt);
            if(isSingleton){
              snvSingCt.addVC(gt);
            }
          }
        subsetList.foreach{ case (sfName,filter,countSetSet) => {
            if(filter.keep(vc)) {
              countSetSet.addVariant(gt=gt,isSingleton=isSingleton,isIndel=isIndel);
            }
        }}
        
      }}
      
      val allTitles : Seq[String] = Seq("HomRef","Het","HomAlt","OtherAlt","Other")
      val altTitles : Seq[String] = Seq("AnyAlt","Het","HomAlt","OtherAlt")
      
      val writer = openWriterSmart(outfile);
      writer.write("sample.ID\t"+"numCalled\tnumMiss\t"+
                   allTitles.mkString("\t")+"\t"+
                   allTitles.map(t => "SNV_" + t).mkString("\t")+"\t"+
                   allTitles.map(t => "INDEL_" + t).mkString("\t")+"\t"+
                   altTitles.map(t => "SING_" + t).mkString("\t")+"\t"+
                   altTitles.map(t => "SINGSNV_" + t).mkString("\t")+"\t"+
                   altTitles.map(t => "SINGINDEL_" + t).mkString("\t")+
                   (subsetList.map{ case (sfName,filter,countSetSet) => { 
                        "\t"+
                        sfName+"_numCalled\t"+sfName+"_numMiss\t"+
                        allTitles.map(t => sfName+"_" + t).mkString("\t")+"\t"+
                        allTitles.map(t => sfName+"_SNV_" + t).mkString("\t")+"\t"+
                        allTitles.map(t => sfName+"_INDEL_" + t).mkString("\t")+"\t"+
                        altTitles.map(t => sfName+"_SING_" + t).mkString("\t")+"\t"+
                        altTitles.map(t => sfName+"_SINGSNV_" + t).mkString("\t")+"\t"+
                        altTitles.map(t => sfName+"_SINGINDEL_" + t).mkString("\t")
                   }}.mkString(""))+
                   "\n")
      Range(0,sampleCt).foreach(i => {
        writer.write(Seq(
            vcfHeader.titleLine.sampleList(i),
            numVariants - varCt.ctMis(i),
            varCt.ctMis(i),
            varCt.getAll(i),
            snvCt.getAll(i),
            indelCt.getAll(i),
            singCt.getAlt(i),
            snvSingCt.getAlt(i),
            indelSingCt.getAlt(i)
        ).mkString("\t") + 
        (subsetList.map{ case (sfName,filter,css) => { 
              "\t"+Seq(
                  css.varCt.getNcalled(i).toString(),
                  css.varCt.ctMis(i),
                  css.varCt.getAll(i),
                  css.snvCt.getAll(i),
                  css.indelCt.getAll(i),
                  css.singCt.getAlt(i),
                  css.snvSingCt.getAlt(i),
                  css.indelSingCt.getAlt(i)
              ).mkString("\t")
        }}.mkString("")) +
        "\n");
      })
      writer.close();
      
      
      
    }
  }

  class RunCalcVariantCountSummary extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "CalcVariantCountSummary", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new BinaryArgument[String](
                                         name = "genotypeTag", 
                                         arg = List("--genotypeTag"), 
                                         valueName = "GT",  
                                         argDesc =  ".",
                                         defaultValue = Some("GT")
                                        ) ::
                    new UnaryArgument( name = "inputFileList",
                                         arg = List("--inputFileList"), // name of value
                                         argDesc = ""+
                                                   "" // description
                                       ) ::
                    new BinaryOptionArgument[Int](
                                         name = "numLinesRead", 
                                         arg = List("--numLinesRead"), 
                                         valueName = "N",  
                                         argDesc =  "Limit file read to the first N lines"
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "subFilterExpressionSets", 
                                         arg = List("--subFilterExpressionSets"), 
                                         valueName = "",  
                                         argDesc =  ""
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "master input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "summaryInfo.txt",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         CalcVariantCountSummary(
             genotypeTag = parser.get[String]("genotypeTag"),
             subFilterExpressionSets = parser.get[Option[String]]("subFilterExpressionSets")
         ).walkVCFFile(
             infile = parser.get[String]("infile"), inputFileList = parser.get[Boolean]("inputFileList"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList"),
             numLinesRead = parser.get[Option[Int]]("numLinesRead")
         )
       }
     }
     
  }
  
  class RedoEnsemblMerge extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "RedoEnsemblMerge", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new UnaryArgument( name = "singleDbFile",
                                         arg = List("--singleDbFile"), // name of value
                                         argDesc = "NOT CURRENTLY SUPPORTED"+
                                                   "" // description
                                       ) ::
                    new BinaryArgument[String](
                                         name = "chromStyle", 
                                         arg = List("--chromStyle"), 
                                         valueName = "hg19",  
                                         argDesc =  ".",
                                         defaultValue = Some("hg19")
                                        ) ::
                    new BinaryOptionArgument[Int](
                                         name = "numLinesRead", 
                                         arg = List("--numLinesRead"), 
                                         valueName = "N",  
                                         argDesc =  "Limit file read to the first N lines"
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "masterCaller", 
                                         arg = List("--masterCaller"), 
                                         valueName = "hc",  
                                         argDesc =  "A caller from which to import ALL info tags."
                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "summaryFile", 
                                         arg = List("--summaryFile"), 
                                         valueName = "summaryData.txt",  
                                         argDesc =  ""
                                        ) :: 
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "variants.vcf",
                                         argDesc = "master input VCF file. Can be gzipped or in plaintext." // description
                                        ) ::
                    new FinalArgument[List[String]](
                                         name = "scVcfFiles",
                                         valueName = "",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[List[String]](
                                         name = "scVcfNames",
                                         valueName = "",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         FixEnsemblMerge2(
             inputVCFs = parser.get[List[String]]("scVcfFiles"),
             inputVcfTypes = parser.get[List[String]]("scVcfNames"),
             masterCaller = parser.get[Option[String]]("masterCaller"),
             summaryFile = parser.get[Option[String]]("summaryFile")
         ).walkVCFFile(
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList"),
             numLinesRead = parser.get[Option[Int]]("numLinesRead")
         )
       }
     }
     
  }
  

  
  case class FixEnsemblMerge2(inputVCFs : Seq[String], inputVcfTypes : Seq[String], masterCaller : Option[String], summaryFile : Option[String]) extends SVcfWalker {

    val LEGAL_VCF_TYPES : Set[String] = Set("hc","fb","ug");
    
    if(inputVcfTypes.exists(! LEGAL_VCF_TYPES.contains(_))){
      error("Illegal VCF TYPE. Must be one of: [\"" + LEGAL_VCF_TYPES.toSeq.sorted.mkString("\",\"") + "\"]");
    }
    LEGAL_VCF_TYPES.foreach(t => {
      if(inputVcfTypes.count(_ == t) > 1){
        error("Illegal VCF TYPES. CANNOT HAVE MORE THAN ONE VCF OF EACH TYPE. FOUND "+inputVcfTypes.count(_ == t) + " VCF's with given type = \""+t+"\"");
      }
    }) 
    
    val fileList = inputVCFs.zip(inputVcfTypes);
    val readers = fileList.map{case (infile,t) => SVcfLine.readVcf(getLinesSmartUnzip(infile), withProgress = false)};
    val headers = readers.map(_._1);
    val iteratorArray : Array[BufferedIterator[SVcfVariantLine]] = readers.map(_._2.buffered).toArray;
    
    val fmtTags = headers.zip(inputVcfTypes).map{ case (h,t) => {
      h.formatLines.map{ fhl => {
           val ct = if(fhl.ID == "AD"){
             "R"
           } else {
             fhl.Number
           }
           (fhl.ID, new SVcfCompoundHeaderLine(in_tag = "FORMAT",ID = t + "_" + fhl.ID, Number = ct, Type = fhl.Type, desc = "For the caller "+t+", " + cleanQuotes(fhl.desc)))
      }}.toSeq
    }}
    
    val genotypeOrdering = new Ordering[String]{
                        def compare(x : String, y : String) : Int = {
                          if(x == y) 0;
                          else if(x == ".") -1;
                          else if(y == ".") 1;
                          else {
                            val xi = string2int(x);
                            val yi = string2int(y);
                            if(xi < yi) -1;
                            else 1;
                          }
                        }
                      }
    
    val masterCallerIdx = masterCaller match {
      case Some(mc) => {
        inputVcfTypes.indexOf(mc)
      }
      case None => None;
    }
    val masterCallerInfoTags : Seq[SVcfCompoundHeaderLine] = masterCaller match {
      case Some(mc) => {
        if( inputVcfTypes.contains(mc)){
          val masterIdx = inputVcfTypes.indexOf(mc)
          val masterHeader = headers(masterIdx);
          masterHeader.infoLines.map{ infoLine => {
            new SVcfCompoundHeaderLine(in_tag = "INFO",ID = "SWH_"+ mc + "_" + infoLine.ID, Number = infoLine.Number, Type = infoLine.Type, desc = "For the caller "+mc+", " + cleanQuotes(infoLine.desc))
          }}
        } else {
          error("ERROR: Master Caller Not Found! Must be one of: " + inputVcfTypes.mkString(","));
          Seq[SVcfCompoundHeaderLine]();
        }
      }
      case None => {
        Seq[SVcfCompoundHeaderLine]();
      }
    }
    
    def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader) = {
      val vcfCodes = VCFAnnoCodes();

      
      val customInfoLines = inputVcfTypes.map{t => {
                                                Seq(
                                                    SVcfCompoundHeaderLine("INFO", vcfCodes.ec_singleCallerAllePrefix+t, ".", "String", "Alt Alleles for caller "+(t))
                                                   )
                                        }}
      

      val extraInfoLines = Seq(
               SVcfCompoundHeaderLine("INFO", vcfCodes.ec_CallMismatch,        "1", "Integer", "Num genotypes that actively disagree (ie called in different ways)"),
               SVcfCompoundHeaderLine("INFO",  vcfCodes.ec_CallMismatchStrict, "1", "Integer", "Num genotypes that do not give the exact same call (including no-call vs call)"),
               SVcfCompoundHeaderLine("INFO",  vcfCodes.ec_EnsembleWarnings,   ".", "String", "List of warnings related to the ensemble calling"),
               SVcfCompoundHeaderLine("INFO",  vcfCodes.ec_alle_callerSets,   "A", "String", "For each alt allele, which callers included the given allele.")
          ) ++ masterCallerInfoTags;
          
          
      val extraFmtLines = Seq(
                SVcfCompoundHeaderLine("FORMAT", "MISMATCH", "1", "String", "All callers do not actively disagree."),
                SVcfCompoundHeaderLine("FORMAT", "MISMATCH_STRICT", "1", "String", "All callers provide the same call."),
                SVcfCompoundHeaderLine("FORMAT", "CallerSupport", ".", "String", "List of callers that support the final call."),
                SVcfCompoundHeaderLine("FORMAT", "HasSupport", "1", "String", "Final call is supported by at least one caller."),
                SVcfCompoundHeaderLine("FORMAT", "ENS_WARN", ".", "String", "Warnings related to the ensemble calls.")
              );
      
      val customFmtLines = inputVcfTypes.map{t =>{
                                                Seq(
                                                    SVcfCompoundHeaderLine("FORMAT", t+"_GT_RAW", ".", "String", "Raw Genotype Call for caller "+t),
                                                    SVcfCompoundHeaderLine("FORMAT", t+"_GT_FIX", "1", "String", "Recoded Genotype Call for caller "+t)
                                                )
                                        }}
      
      val newFmtLines = vcfHeader.formatLines ++ fmtTags.flatMap{ newTagLines =>
                                        newTagLines.map(_._2)
      } ++ customFmtLines.flatten ++ extraFmtLines;
      
      val newHeader = SVcfHeader(infoLines = vcfHeader.infoLines ++ customInfoLines.flatten ++ extraInfoLines, 
                                formatLines = newFmtLines,
                                otherHeaderLines = vcfHeader.otherHeaderLines,
                                titleLine = vcfHeader.titleLine);
      
      val sampNames = vcfHeader.titleLine.sampleList;
      val sampCt = sampNames.length;
      //var currIter = vcIter.buffered;
      
      /*
      val out = groupBySpan(vcIter)(vc => vc.pos).flatMap( vcSeq => {
        val currPos = vcSeq.head.pos;
        iteratorArray.indices.foreach{i => {
          //iteratorArray(i) = iteratorArray(i).dropWhile(vAlt => vAlt.pos < currPos);
          
        }}
        val otherVcAtPos = iteratorArray.indices.map{i => {
          val (a,r) = spanVector(iteratorArray(i))(vAlt => vAlt.pos == currPos);
          iteratorArray(i) = r;
          a
        }}*/
      
      val out = groupBySpan(vcIter.buffered)(vc => vc.pos).flatMap( vcSeq => {
        val currPos = vcSeq.head.pos;
        val currChrom = vcSeq.head.chrom;
        iteratorArray.indices.foreach{i => {
          //iteratorArray(i) = iteratorArray(i).dropWhile(vAlt => vAlt.pos < currPos);
          skipWhile(iteratorArray(i))(vAlt => vAlt.chrom != currChrom);
          skipWhile(iteratorArray(i))(vAlt => vAlt.pos < currPos && vAlt.chrom == currChrom);
        }}
        val otherVcAtPos = iteratorArray.indices.map{i => {
          extractWhile(iteratorArray(i))(vAlt => vAlt.pos == currPos && vAlt.chrom == currChrom);
        }}
        vcSeq.iterator.map(vc => {
          
          val vb = vc.getOutputLine();
          val altAlles = vc.alt;
          var ensembleWarnings = Set[String]();
          val ploidy = vc.genotypes.genotypeValues(0).map{ _.split("/").length }.max
          
          if(ploidy > 2){
            warning("Ploidy greater than 2! Ploidy = "+ploidy,"POLYPLOID",100)
          } else if(ploidy == 1){
            warning("Haploid! Ploidy = "+ploidy,"HAPLOID",100)
          }
          
          if(altAlles.length > 1){
            ensembleWarnings = ensembleWarnings + ("MULTIALLELIC");
          }
          
          var callerSets = Array.fill[Set[String]](altAlles.length)(Set[String]());
          var sampleWarn = Array.fill[Set[String]](sampCt)(Set[String]()); 
          
          try{
          otherVcAtPos.zip(fmtTags).zipWithIndex.foreach{ case ((otherLines,otherFmtTags),otherFileIdx) => {
            val otherFileType = inputVcfTypes(otherFileIdx);
            val matchIdx = otherLines.zipWithIndex.flatMap{ case (otherVC, otherLineIdx) => {
              otherVC.alt.zipWithIndex.filter{case (a,idx) => {a != "*"}}.flatMap{ case (otherAlle, otherAlleIdx) => {
                vb.alt.zipWithIndex.filter{case (alle,idx) => {alle == otherAlle}}.map{case (currAlle,currAlleIdx) => ((currAlleIdx,(otherLineIdx,otherAlleIdx)))}
              }}
            }}.toMap;
            val linesWithMatch = matchIdx.map{case (currAlleIdx,(otherLineIdx,otherAlleIdx)) => otherLineIdx}.toSet.toSeq.sorted;
            val linesWithoutMatch = (otherLines.indices.toSet -- linesWithMatch.toSet).toVector.sorted;
            val numLinesWithMatch = linesWithMatch.size;
            
            if(otherLines.length > 1){
              warning("Multiple lines found at location (Caller "+inputVcfTypes(otherFileIdx)+", POS="+vcSeq.head.chrom+":"+currPos+")","MULTILINE_LOCUS_"+otherFileType,100);
              ensembleWarnings = ensembleWarnings + ("MULTILINELOCUS_"+inputVcfTypes(otherFileIdx));
            }
            
            if(numLinesWithMatch > 1){
              warning("Multiple lines found at location that contain matches (Caller "+inputVcfTypes(otherFileIdx)+", POS="+vcSeq.head.chrom+":"+currPos+")","MULTIMATCHLINE_LOCUS_"+otherFileType,100);
            }
            
            if(! matchIdx.isEmpty){
              val alts = altAlles.zipWithIndex.filter{case (a,idx) => {a != "*"}};
              val fmtA = fmtTags(otherFileIdx).filter{case (rawFmtTag,fmtLine) => {
                fmtLine.Number == "A"
              }}.map{ case (rawFmtTag,fmtLine) => {
                (Array.fill[String](sampCt,alts.length)("."),rawFmtTag,fmtLine);
              }}
              val fmtR = fmtTags(otherFileIdx).filter{case (rawFmtTag,fmtLine) => {
                fmtLine.Number == "R"
              }}.map{ case (rawFmtTag,fmtLine) => {
                (Array.fill[String](sampCt,alts.length + 1)("."),rawFmtTag,fmtLine);
              }}
              val fmtOther = fmtTags(otherFileIdx).filter{case (rawFmtTag,fmtLine) => {
                fmtLine.Number != "R" && fmtLine.Number != "A";
              }}.map{ case (rawFmtTag,fmtLine) => {
                (Array.fill[String](sampCt,numLinesWithMatch)("."),rawFmtTag,fmtLine);
              }}
              
              val altAlleArray = Array.fill[String](otherLines.length)(".");
              val altGtArray = Array.fill[String](sampCt,otherLines.length)("./.");
              val altGtFixedArray = Array.fill[String](sampCt,2)(".");
              
              matchIdx.groupBy{ case (currAlleIdx,(otherLineIdx,otherAlleIdx)) => { otherLineIdx }}.foreach{ case (otherLineIdx, s) => {
                
                val matchLineIdx = linesWithMatch.indexOf(otherLineIdx);
                val otherAlts = otherLines(otherLineIdx).alt.zipWithIndex
                try{
                  val otherGenotypeFormat = otherLines(otherLineIdx).format;
                  val otherGenotypeArray = otherLines(otherLineIdx).genotypes.genotypeValues;
                  val rawArrayFmtA = fmtA.zipWithIndex.filter{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    otherGenotypeFormat.contains(rawFmtTag);
                  }}.map{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    val fmtIdx = otherGenotypeFormat.indexOf(rawFmtTag);
                    (otherGenotypeArray(fmtIdx).map(otherG => {
                      if(otherG == ".") {
                        Array.fill[String](otherAlts.length)(".");
                      } else {
                        otherG.split(",")
                      }
                    }),i);
                  }}
                  val rawArrayFmtR = fmtR.zipWithIndex.filter{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    otherGenotypeFormat.contains(rawFmtTag);
                  }}.map{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    val fmtIdx = otherGenotypeFormat.indexOf(rawFmtTag);
                    (otherGenotypeArray(fmtIdx).map(otherG => {
                      if(otherG == ".") {
                        Array.fill[String](otherAlts.length + 1)(".");
                      } else {
                        otherG.split(",")
                      }
                    }),i);
                  }}
                  val rawArrayFmtOther = fmtOther.zipWithIndex.filter{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    otherGenotypeFormat.contains(rawFmtTag);
                  }}.map{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    val fmtIdx = otherGenotypeFormat.indexOf(rawFmtTag);
                    (otherGenotypeArray(fmtIdx),i);
                  }}
                  
                  fmtOther.zipWithIndex.filter{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    otherGenotypeFormat.contains(rawFmtTag);
                  }}.foreach{ case ((gArray,rawFmtTag,fmtLine),i) => {
                    val fmtIdx = otherGenotypeFormat.indexOf(rawFmtTag);
                    val gArray = otherGenotypeArray(fmtIdx)
                    Range(0,sampCt).foreach{sampIdx => {
                      val currArray = fmtOther(i)._1;
                      currArray(sampIdx)(matchLineIdx) = gArray(sampIdx);
                    }}
                  }}
                  
                  rawArrayFmtR.foreach{ case (gArray,i) => {
                      val currArray = fmtR(i)._1;
                      Range(0,sampCt).foreach{sampIdx => {
                        currArray(sampIdx)(0) = gArray(sampIdx)(0);
                      }}
                  }}
                  
                  s.foreach{ case (currAlleIdx,(oli,otherAlleIdx)) => {
                    callerSets(currAlleIdx) = callerSets(currAlleIdx) + inputVcfTypes(otherFileIdx);
                    rawArrayFmtA.foreach{ case (gArray,i) => {
                      val currArray = fmtA(i)._1;
                      Range(0,sampCt).foreach{sampIdx => {
                        if( otherAlleIdx >= gArray(sampIdx).length) warning("Found problem with genotype "+sampIdx+": otherAlleIdx="+otherAlleIdx+", but for tag "+fmtA(i)._2+", length is: "+gArray(sampIdx).length,"MALFORMED_GENOTYPE",100);
                        currArray(sampIdx)(currAlleIdx) = gArray(sampIdx)(otherAlleIdx);
                      }}
                    }}
                    rawArrayFmtR.foreach{ case (gArray,i) => {
                      val currArray = fmtR(i)._1;
                      Range(0,sampCt).foreach{sampIdx => {
                        currArray(sampIdx)(currAlleIdx + 1) = gArray(sampIdx)(otherAlleIdx + 1);
                      }}
                    }}
                    Range(0,sampCt).foreach{sampIdx => {
                      val geno = otherGenotypeArray(0)(sampIdx).split("/");
                      geno.zipWithIndex.foreach{case (g,i) => {
                        if(g == (otherAlleIdx + 1).toString()){
                          if(altGtFixedArray(sampIdx)(i) != "." && altGtFixedArray(sampIdx)(i) != (currAlleIdx + 1).toString){
                              warning("Overwriting existing variant on a multiline merger!","OVERWRITE_VARIANT_ON_MULTILINE_MERGE_"+otherFileType,100);
                              ensembleWarnings = ensembleWarnings + ("OVERWRITE_VARIANT_ON_MULTILINE_MERGE_"+otherFileType);
                              sampleWarn(sampIdx) = sampleWarn(sampIdx) + ("OVERWRITE_VARIANT_ON_MULTILINE_MERGE_"+otherFileType);
                          }
                          altGtFixedArray(sampIdx)(i) = (currAlleIdx + 1).toString;
                        }
                      }}
                    }}
                  }}
                  altAlleArray(otherLineIdx) = otherLines(otherLineIdx).alt.mkString(",");
                  Range(0,sampCt).foreach{ sampIdx => {
                    altGtArray(sampIdx)(otherLineIdx) = otherGenotypeArray(0)(sampIdx);
                  }}
                  Range(0,sampCt).foreach{sampIdx => {
                      val geno = otherGenotypeArray(0)(sampIdx).split("/");
                      geno.zipWithIndex.foreach{case (g,i) => {
                        if(g == "0"){
                          if(altGtFixedArray(sampIdx)(i) != "." && altGtFixedArray(sampIdx)(i) != "0"){
                              warning("Overwriting existing variant on a multiline merger!","OVERWRITE_REFVARIANT_ON_MULTILINE_MERGE_"+otherFileType,100);
                              ensembleWarnings = ensembleWarnings + ("OVERWRITE_REFVARIANT_ON_MULTILINE_MERGE_"+otherFileType);
                              sampleWarn(sampIdx) = sampleWarn(sampIdx) + ("OVERWRITE_VARIANT_ON_MULTILINE_MERGE_"+otherFileType);
                          } else {
                            altGtFixedArray(sampIdx)(i) = "0";
                          }
                        }
                      }}
                  }}
                } catch {
                  case e : Exception => {
                    reportln("Caught exception in matchIDX iteration.\n"+
                       "###CURRENT OTHER VCF LINE:\n"+
                       otherLines(otherLineIdx).getVcfString + "\n"+
                       s.map{ case (currAlleIdx,(oli,otherAlleIdx)) => {
                          "   (currAlleIdx="+currAlleIdx+", otherLineIdx="+otherLineIdx+", otherAlleIdx="+otherAlleIdx+", matchLineIdx="+matchLineIdx+")"
                       }}.mkString("\n")+"\n"+
                       
                       "###Match IDX:"+
                       matchIdx.map{ case (cai,(oli,oai)) => {
                         "   (currAlleIdx="+cai+", otherLineIdx="+oli+", otherAlleIdx="+oai+")\n"+
                         otherLines(oli).getVcfString
                       }}.mkString("\n")+"\n"+
                       "","note");
                    throw e;
                  }
                }
              }}
              
              linesWithoutMatch.foreach{otherLineIdx => {
                val otherGenotypeArray = otherLines(otherLineIdx).genotypes.genotypeValues;
                Range(0,sampCt).foreach{sampIdx => {
                    val geno = otherGenotypeArray(0)(sampIdx).split("/");
                    geno.zipWithIndex.foreach{case (g,i) => {
                      if(g == "0"){
                        if(altGtFixedArray(sampIdx)(i) != "." && altGtFixedArray(sampIdx)(i) != "0"){
                            warning("Conflicting existing refvariant on a multiline merger!","OVERWRITE_REFVARIANT_ON_MULTILINE_MERGE_"+otherFileType,100);
                            ensembleWarnings = ensembleWarnings + ("OVERWRITE_REFVARIANT_ON_MULTILINE_MERGE_"+otherFileType);
                            sampleWarn(sampIdx) = sampleWarn(sampIdx) + ("OVERWRITE_VARIANT_ON_MULTILINE_MERGE_"+otherFileType);
                        } else {
                          altGtFixedArray(sampIdx)(i) = "0";
                        }
                      }
                    }}
                }}
              }}
              
              if(ploidy > 1){
                Range(0,sampCt).foreach{sampIdx => {
                      altGtFixedArray(sampIdx) = altGtFixedArray(sampIdx).sortBy(s => s)(genotypeOrdering)
                }}
              }
              //vb.genotypes.fmt = vb.genotypes.fmt ++ fmtA.map{_._3} ++ fmtR.map{_._3} ++ fmtOther.map{_._3} ++ customFmtLines(otherFileIdx)
              vb.genotypes.fmt = vb.genotypes.fmt ++ fmtA.map{_._3.ID} ++ fmtR.map{_._3.ID} ++ fmtOther.map{_._3.ID} ++ customFmtLines(otherFileIdx).map(_.ID);
              //vb.in_format = vb.in_format ++ fmtA.map{_._3.ID} ++ fmtR.map{_._3.ID} ++ fmtOther.map{_._3.ID} ++ customFmtLines(otherFileIdx).map(_.ID);
              vb.genotypes.genotypeValues = ( vb.genotypes.genotypeValues ++ fmtA.map{_._1.map(_.mkString(","))} ++ 
                                              fmtR.map{_._1.map(_.mkString(","))} ++ 
                                              fmtOther.map{_._1.map(_.mkString("|"))} ) ++
                                             Array(
                                                 altGtArray.map(_.mkString("|")),
                                                 altGtFixedArray.map(_.mkString("/"))
                                             )
              vb.in_info = vb.in_info ++ Map(
                    (vcfCodes.ec_singleCallerAllePrefix+otherFileType,Some(altAlleArray.mkString("|")))
                  )
              if(masterCaller.isDefined && masterCaller.get == inputVcfTypes(otherFileIdx)){
                 if(linesWithMatch.size > 1){
                   warning("WARNING: Master Caller should only ever have 1 line per position. Ambiguous output!","BAD_MASTER_CALLER",100);
                 }
                 val otherLineIdx = linesWithMatch.head;
                 val otherInfo = otherLines(otherLineIdx).info;
                 vb.in_info = vb.in_info ++ otherInfo.map{ case (infoTag,infoValue) => {
                   ("SWH_"+ masterCaller.get + "_" + infoTag , infoValue)
                 }}
              }
            }  
          }}
          
          
          val masterGT = vb.genotypes.genotypeValues(0);
          val mm = Array.fill[Boolean](sampCt)(false);
          val mmStrict = Array.fill[Boolean](sampCt)(false);
          inputVcfTypes.foreach{ivt => {
              val gtTag = ivt+"_GT_FIX";
              val gtIdx = vb.format.indexOf(gtTag);
              if(gtIdx != -1){
                val otherGT = vb.genotypes.genotypeValues(gtIdx);
                masterGT.indices.foreach{sampIdx => {
                  mmStrict(sampIdx) = mmStrict(sampIdx) || masterGT(sampIdx) != otherGT(sampIdx);
                  val mgt = masterGT(sampIdx).split("/");
                  val ogt = otherGT(sampIdx).split("/");
                  mgt.zip(ogt).foreach{case (m,o) => {
                    if(m != "." && o != "." && m != o){
                      mm(sampIdx) = true;
                      sampleWarn(sampIdx) = sampleWarn(sampIdx) + ("CALLER_MISMATCH_"+ivt);
                    }
                  }}
                }}
              }
          }}
          val gtTags = inputVcfTypes.map{ivt => { ivt+"_GT_FIX" }};
          val gtSet  = gtTags.zip(gtTags.map{vb.format.indexOf(_)}).zip(inputVcfTypes).filter{ case ((t,i),f) => i != -1};
          
          val gtCallerSupport = masterGT.indices.map{sampIdx => {
            val mgt = masterGT(sampIdx);
            if(! mgt.contains('.')){
              gtSet.flatMap{ case ((otherTag,otherIdx),ivt) => {
                if(vb.genotypes.genotypeValues(otherIdx)(sampIdx) == mgt){
                  Some(ivt);
                } else {
                  None;
                }
              }}.mkString(",");
            } else {
              "NA"
            }
          }}.toArray;
          val gtIsSupported = masterGT.indices.map{sampIdx => {
            val mgt = masterGT(sampIdx).split("/");
            val isSupported = mgt.forall{m => {
              m == "." || {
                gtSet.exists{ case ((otherTag,otherIdx),ivt) => {
                  val ogt = vb.genotypes.genotypeValues(otherIdx)(sampIdx).split("/");
                  ogt.contains(m);
                }}
              }
            }}
            if(! isSupported){
              sampleWarn(sampIdx) = sampleWarn(sampIdx) + ("UNSUPPORTED");
            }
            isSupported
          }}
          val numUnsupportedGT = gtIsSupported.count(! _);
          if(numUnsupportedGT > 0){
            val warnMsg = "UNSUPPORTED_GT"
            warning("Found "+numUnsupportedGT+ " unsupported genotypes!",warnMsg,100);
            ensembleWarnings = ensembleWarnings + (warnMsg);
          }
          
          //vb.genotypes.fmt = vb.genotypes.fmt ++ extraFmtLines
          vb.genotypes.fmt = vb.genotypes.fmt ++ extraFmtLines.map{efl => efl.ID}
          //vb.in_format = vb.in_format ++ extraFmtLines.map{efl => efl.ID}
          vb.genotypes.genotypeValues = vb.genotypes.genotypeValues ++ 
                                               Array(mm.map(if(_) "1" else "0"),
                                                     mmStrict.map(if(_) "1" else "0"),
                                                     gtCallerSupport,
                                                     gtIsSupported.map(if(_) "1" else "0").toArray,
                                                     sampleWarn.map{s => s.toSeq.padTo(1,".").sorted.mkString(",")})
          
                           
          vb.in_info = vb.in_info ++ Map(
                ( vcfCodes.ec_CallMismatch,Some( mm.count(x => x).toString )),
                ( vcfCodes.ec_CallMismatchStrict,Some( mmStrict.count(x => x).toString )),
                (vcfCodes.ec_EnsembleWarnings,Some( ensembleWarnings.toSeq.sorted.padTo(1,".").mkString(",") )),
                (vcfCodes.ec_alle_callerSets,Some(callerSets.map{s => s.toSeq.padTo(1,".").sorted.mkString("|")}.padTo(1,".").mkString(",")))
              );
          
          } catch {
            case e : Exception => {
              reportln("Caught exception in VCF iteration.\n"+
                       "###MASTER VCF LINE:\n"+
                       vc.getVcfString+"\n"+
                       (if(vcSeq.length > 1){ "###All master lines at position:\n" + vcSeq.map(v => v.getVcfString).mkString("\n") + "\n" } else { "" }) +
                       "###All VCF Lines at position:\n"+
                       otherVcAtPos.zip(inputVcfTypes).map{ case (otherVCs,fileID) => {
                         otherVCs.map(v => fileID + "\t" + v.getVcfString).mkString("\n")
                       }}.mkString("\n") + "\n"+
                       "","note");
              
              throw e;
            }
          }
          vb;
          
        })
      })
      
      (out,newHeader);
    };
    
  }
  
  
  class CommandFilterVCF extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "FilterVCF", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new BinaryOptionArgument[Int](
                                         name = "numLinesRead", 
                                         arg = List("--numLinesRead"), 
                                         valueName = "N",  
                                         argDesc =  "Limit file read to the first N lines"
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile.vcf.gz",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "variantFilterExpression",
                                         valueName = "expr",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file or a comma-delimited list of files. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         VcfExpressionFilter(
             filterExpr = parser.get[String]("variantFilterExpression")
         ).walkVCFFile(
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList"),
             numLinesRead = parser.get[Option[Int]]("numLinesRead")
         )
       }
     }
  }

  class CommandVcfToMatrix extends CommandLineRunUtil {
     override def priority = 20;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "VcfToMatrix", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "" + ALPHA_WARNING,
          argList = 
                    new BinaryOptionArgument[Int](
                                         name = "numLinesRead", 
                                         arg = List("--numLinesRead"), 
                                         valueName = "N",  
                                         argDesc =  "Limit file read to the first N lines"
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "chromList", 
                                         arg = List("--chromList"), 
                                         valueName = "chr1,chr2,...",  
                                         argDesc =  "List of chromosomes. If supplied, then all analysis will be restricted to these chromosomes. All other chromosomes wil be ignored."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "variantFile", 
                                         arg = List("--variantFile"), 
                                         valueName = "varfile.txt",  
                                         argDesc =  ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "variantFilterExpression", 
                                         arg = List("--variantFilterExpression"), 
                                         valueName = "...",  
                                         argDesc =  "Variant-level filter expression."
                                        ) ::
                    new UnaryArgument( name = "variantsAsRows",
                                         arg = List("--variantsAsRows"), // name of value
                                         argDesc = "..."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "variantsAsColumns",
                                         arg = List("--variantsAsColumns"), // name of value
                                         argDesc = "Default behavior, this parameter has no effect."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "writeTitleColumn",
                                         arg = List("--writeTitleColumn"), // name of value
                                         argDesc = "..."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "writeHeaderRow",
                                         arg = List("--writeHeaderRow"), // name of value
                                         argDesc = "..."+
                                                   "" // description
                                       ) ::
                    new BinaryArgument[String](name = "naString",
                                           arg = List("--naString"),  
                                           valueName = "NA", 
                                           argDesc = "", 
                                           defaultValue = Some("NA")
                                           ) :: 
                    new BinaryArgument[String](name = "otherString",
                                           arg = List("--otherString"),  
                                           valueName = "NA", 
                                           argDesc = "", 
                                           defaultValue = Some("NA")
                                           ) :: 
                    new BinaryArgument[String](name = "delimString",
                                           arg = List("--delimString"),  
                                           valueName = "\\t", 
                                           argDesc = "", 
                                           defaultValue = Some("\t")
                                           ) :: 
                    new BinaryArgument[String](name = "gtTagString",
                                           arg = List("--gtTagString"),  
                                           valueName = "\\t", 
                                           argDesc = "", 
                                           defaultValue = Some("GT")
                                           ) ::    
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile.vcf.gz",
                                         argDesc = "" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.vcf.gz",
                                         argDesc = "The output file. Can be gzipped or in plaintext."// description
                                        ) ::
                    internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );

     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       if(out){
         VcfToMatrix(
             variantsAsRows = parser.get[Boolean]("variantsAsRows"),
             columnDelim = parser.get[String]("delimString"),
             naString = parser.get[String]("naString"),
             multiAlleString = parser.get[String]("otherString"),
             writeHeader = parser.get[Boolean]("writeHeaderRow"),
             writeTitleColumn = parser.get[Boolean]("writeTitleColumn"),
             variantFile = parser.get[Option[String]]("variantFile"),
             gtTagString = parser.get[String]("gtTagString")
         ).walkVCFFile(
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             chromList = parser.get[Option[List[String]]]("chromList"),
             numLinesRead = parser.get[Option[Int]]("numLinesRead"),
             variantFilterExpression = parser.get[Option[String]]("variantFilterExpression")
         )
       }
     }
  }
  
  
  case class VcfToMatrix(variantsAsRows : Boolean = false, columnDelim : String = "\t", 
                         naString : String = "NA", multiAlleString : String = "NA",
                         writeHeader : Boolean = false, writeTitleColumn : Boolean = false,
                         variantFile : Option[String] = None,
                         gtTagString : String = "GT")  {
    
    def walkVCFFile(infile : String, outfile : String, 
                    variantFilterExpression : Option[String] = None,
                    chromList : Option[List[String]],numLinesRead : Option[Int],
                    allowVcfList : Boolean = true){
      
      val indata = if(infile.contains(',') && allowVcfList){
        val infiles = infile.split(",");
        val allInputLines = flattenIterators(infiles.iterator.map{inf => addIteratorCloseAction(iter =getLinesSmartUnzip(inf), closeAction = (() => {reportln("finished reading file: "+inf,"note")}))}).buffered
        val headerLines = extractWhile(allInputLines)( a => a.startsWith("#"));
        val remainderLines = allInputLines.filter( a => ! a.startsWith("#"));
        headerLines.iterator ++ remainderLines;
      } else {
        getLinesSmartUnzip(infile)
      }
      
      val (vcfHeader,vcIter) = if(chromList.isEmpty){
        SVcfLine.readVcf(indata,withProgress = true)
      } else if(chromList.get.length == 1){
        val chrom = chromList.get.head;
        SVcfLine.readVcf(indata.filter{line => {
          line.startsWith(chrom+"\t") || line.startsWith("#")
        }},withProgress = true)
      } else {
        val chromSet = chromList.get.toSet;
        val (vh,vi) = SVcfLine.readVcf(indata,withProgress = true)
        (vh,vi.filter(line => { chromSet.contains(line.chrom) }))
      } 
      
      
      val (vcIter2,vcfHeader2) = variantFilterExpression match {
        case Some(expr) => {
          VcfExpressionFilter(expr).walkVCF(vcIter,vcfHeader)
        }
        case None => {
          (vcIter,vcfHeader)
        }
      }
      
      val vcIter3 = if(numLinesRead.isDefined){
        vcIter2.take(numLinesRead.get);
      } else {
        vcIter2
      }
      
      val lineIter = walkVCF(vcIter3,vcfHeader2,verbose=true);
      
      var lnct = 0;
      val writer = openWriterSmart(outfile);
      lineIter.foreach{line => {
        lnct = lnct + 1;
        writer.write(line+"\n");
      }}
      writer.close();
      reportln("Wrote "+lnct + " lines.","note");
    }
    
    def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : Iterator[String] = {
      val sampleList = vcfHeader.titleLine.sampleList
      
      if(variantsAsRows){
        (
          if(writeHeader){
            Iterator[String](sampleList.mkString(columnDelim))
          } else {
            Iterator[String]();
          }
        ) ++ vcIter.zipWithIndex.map{ case (vc,i) => {
          val gtIdx = vc.format.indexOf(gtTagString);
          (if(writeTitleColumn){ i+columnDelim } else {""}) + vc.genotypes.genotypeValues(gtIdx).map{ g => {
            if(g == "./." || g == "."){
              "."
            } else if(g == "0/0"){
              "0"
            } else if(g == "0/1"){
              "1"
            } else if(g == "1/1"){
              "2"
            } else {
              naString
            }
          }}.mkString(columnDelim);
        }}
      } else {
        val (genoArray,varArray) = vcIter.map{ vc => {
          val gtIdx = vc.format.indexOf(gtTagString);
          (vc.genotypes.genotypeValues(gtIdx).map{ g => {
            if(g == "./." || g == "."){
              naString
            } else if(g == "0/0" || g == "0"){
              "0"
            } else if(g == "0/1"){
              "1"
            } else if(g == "1/1" || g == "1"){
              "2"
            } else {
              multiAlleString
            }
          }},
            vc.getSimpleVcfString()
          )
        }}.toVector.unzip
        if(genoArray.isEmpty){
          warning("Writing 0 variants!","WRITING_ZERO_VARIANTS",-1);
          Iterator[String]();
        } else {
          reportln("Writing "+genoArray.length+" variants, "+genoArray(0).length+" samples.","note");
          variantFile match {
            case Some(f) => {
              val writer = openWriterSmart(f);
              varArray.foreach{v => {
                writer.write(v + "\n");
              }}
              writer.close();
            }
            case None => {
              //do nothing
            }
          }
          (
            if(writeHeader){
              Iterator[String](genoArray(0).indices.mkString(columnDelim))
            } else {
              Iterator[String]();
            }
          ) ++ genoArray(0).indices.iterator.map{ j => {
            (if(writeTitleColumn){ sampleList(j) + columnDelim} else {""})+ genoArray.indices.map{ i => {
              genoArray(i)(j)
            }}.mkString(columnDelim)
          }}
        }
        
      }
    }
  }
  
  case class VcfExpressionFilter(filterExpr : String, explainFile : Option[String] = None) extends SVcfWalker {

    /*def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader) = {
       (vcIter, vcfHeader)
    }*/
    
    val parser : SVcfFilterLogicParser = internalUtils.VcfTool.SVcfFilterLogicParser();
    val filter : SFilterLogic[SVcfVariantLine] = parser.parseString(filterExpr);
    //val parser = internalUtils.VcfTool.SVcfFilterLogicParser();
    //val filter = parser.parseString(filterExpr);
    reportln("Parsed filter:\n   "+filter.printTree(),"note");
    
    def walkVCF(vcIter : Iterator[SVcfVariantLine], vcfHeader : SVcfHeader, verbose : Boolean = true) : (Iterator[SVcfVariantLine],SVcfHeader) = {
       (vcIter.filter{ vc => {
         filter.keep(vc);
       }}, vcfHeader)
    }
  }
  
}




























