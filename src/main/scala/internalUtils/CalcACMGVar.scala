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
     override def priority = 30;
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
                                         name = "canonicalTxFile", 
                                         arg = List("--canonicalTxFile"), 
                                         valueName = "knownCanonical.txt",
                                         //argDesc =  "A file containing a list of transcript ID's and whether or not they are a RefSeq transcript. Must have at least 2 labelled columns, name and isRefSeq. The header line may begin with a #, or not. It can be compressed or in plaintext."
                                         argDesc =  "A file containing a list of canonical transcript IDs. It must have a header line with a column labelled \"transcript\". The header line may begin with a #, or not. It can be compressed or in plaintext."
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
                                         name = "inSilicoKeys", 
                                         arg = List("--inSilicoKeys"), 
                                         valueName = "dbNSFP_MetaSVM_pred:D:T",
                                         argDesc =  "This must be a comma-delimited list (no spaces). Each element in the list must consist of 3 parts, seperated by colons (\":\"). "+
                                                    "The first part is the INFO key referring to the stored results of an in silico prediction algorithm. The 2nd is a \"|\"-delimited list "+
                                                    "of the values that should be interpreted as \"predicted damaging\", and the third is a \"|\"-delimited list of the values that should be interpreted as \"predicted benign\""
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
                    new BinaryOptionArgument[String](
                                         name = "summaryOutputFile", 
                                         arg = List("--summaryOutputFile"), 
                                         valueName = "summaryOutputFile.txt",  
                                         argDesc =  "Optional summary output file."
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
                   rmskFile = parser.get[Option[String]]("rmskFile"),
                   refSeqFile = parser.get[Option[String]]("canonicalTxFile"),
                   inSilicoParams = parser.get[Option[String]]("inSilicoKeys")
                   //dropKeys = parser.get[Option[List[String]]]("dropKeys"),
                   ).chain(CalcACMGVar.SummaryACMGWalker(
                             groupFile = parser.get[Option[String]]("groupFile"),
                             groupList = None,
                             superGroupList = parser.get[Option[String]]("superGroupList"),
                             outfile = parser.get[Option[String]]("summaryOutputFile").getOrElse("undefinedfile.txt")
                           ), flag = parser.get[Option[String]]("summaryOutputFile").isDefined
                   ).walkVCFFile(
                     infile    = parser.get[String]("invcf"),
                     outfile   = parser.get[String]("outvcf"),
                     chromList = parser.get[Option[List[String]]]("chromList")
                   )
           
     }
  }
  }
  
  
  
  /*
   * 
         Vector[(String,Set[String],Set[String])](
          ("dbNSFP_MetaSVM_pred",Set[String]("D"),Set[String]("T"))//, //Includes predictors:   PolyPhen-2, SIFT, LRT, MutationTaster, Mutation Assessor, FATHMM, GERP++, PhyloP, SiPhy
                        //("dbNSFP_PROVEAN_pred",Set[String]("D"))
                    )
   * 
   */
  case class SummaryACMGWalker(
                   groupFile : Option[String],
                   groupList : Option[String],
                   superGroupList : Option[String],
                   outfile : String
                 ) extends internalUtils.VcfTool.VCFWalker {
      var vcfCodes : VCFAnnoCodes = VCFAnnoCodes();
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
      val groups = groupSet.toVector.sorted;
      reportln("Final Groups:","debug");
      for((g,i) <- groups.zipWithIndex){
          reportln("Group "+g+" ("+groupToSampleMap(g).size+")","debug");
      }
      
      val geneAndTypeCountMap = new scala.collection.mutable.AnyRefMap[String,Map[String,Int]](x => Map[String,Int]().withDefault(s => 0));
      val geneAndTypeCountMap_CANON = new scala.collection.mutable.AnyRefMap[String,Map[String,Int]](x => Map[String,Int]().withDefault(s => 0));
      
      
     /* val variantSampleMap = new scala.collection.mutable.AnyRefMap[Int,Set[String]](x => Set[String]());
      val variantGeneMap = new scala.collection.mutable.AnyRefMap[Int,Set[String]]();
      val variantTypeMap = new scala.collection.mutable.AnyRefMap[Int,String]();
      val variantTypeMap_CANON = new scala.collection.mutable.AnyRefMap[Int,String]();*/
      
      //(geneSet, (level,level_canon), sampleSet)
      var rawVariantInfo = Vector[(Set[String],(String,String),Set[String])]();
      var allGeneSet = Set[String]();
      
      def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
        val sampNames = vcfHeader.getSampleNamesInOrder().toVector;
        
        val newHeaderLines = List[VCFHeaderLine](); 
        val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
        
        return (addIteratorCloseAction(vcIter.zipWithIndex.map{ case (vc,i) => {
           val geneList = getAttributeSmart(vc,vcfCodes.assess_geneList).toSet.toVector.toSet;
           val lvl = vc.getAttributeAsString(vcfCodes.assess_RATING,"UNK");
           val lvl_CANON = vc.getAttributeAsString(vcfCodes.assess_RATING_CANON,"UNK");
           
           allGeneSet = allGeneSet ++ geneList;
           
           val refAlle = vc.getReference();
           val altAlleles = Range(0,vc.getNAlleles()-1).map((a) => vc.getAlternateAllele(a)).zipWithIndex.filter{case (alt,altIdx) => { alt.getBaseString() != "*" }}
           if( altAlleles.length != 1 ){
             error("More than one alt allele found!");
           }
           if( altAlleles.head._2 != 0){
             warning("Alt alleIdx != 0","ALTIDX_NE_0",20);
           }
           val alt = altAlleles.head._1;
           
           //variantGeneMap(i) = geneList.toSet;
           //variantTypeMap(i) = 
           val altGenoSampleSet = vc.getGenotypes().iterator().filter(g => {
             g.getAlleles().exists(a => a.equals(alt))
           }).map(g => {
             g.getSampleName();
           }).toSet;
           
           if(! vc.isFiltered()){
             rawVariantInfo = rawVariantInfo :+ (geneList,(lvl,lvl_CANON), altGenoSampleSet);
           }
           
           vc;
        }},closeAction=() => {
          
          reportln("Filtering VariantInfo..." + getDateAndTimeString,"debug");
          val variantInfo = rawVariantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
            lvl == "PATHO" || lvl == "LPATH" || lvl_CANON == "PATHO" || lvl_CANON == "LPATH";
          }}
          reportln("Done filtering VariantInfo. " + getDateAndTimeString,"debug");
          val out = openWriterSmart(outfile);
          val allGenes = allGeneSet.toVector.sorted;
          
          reportln("   geneSampCount_patho..." + getDateAndTimeString,"debug");
          val geneSampCount_patho = allGenes.map{ gene => { 
            val geneVars = variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              geneSet.contains(gene) && lvl == "PATHO"
            }}
            groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            geneVars.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}}}
          reportln("   geneSampCount_patho_canon..." + getDateAndTimeString,"debug");
          val geneSampCount_patho_canon = allGenes.map{ gene => { 
            val geneVars = variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              geneSet.contains(gene) && lvl_CANON == "PATHO"
            }}
            groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            geneVars.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}}}
          reportln("   geneSampCount_lpath..." + getDateAndTimeString,"debug");
          val geneSampCount_lpath = allGenes.map{ gene => { 
            val geneVars = variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              geneSet.contains(gene) && (lvl == "PATHO" || lvl == "LPATH")
            }}
            groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            geneVars.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}}}
          reportln("   geneSampCount_lpath_canon..." + getDateAndTimeString,"debug");
          val geneSampCount_lpath_canon = allGenes.map{ gene => { 
            val geneVars = variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              geneSet.contains(gene) && (lvl_CANON == "PATHO" || lvl_CANON == "LPATH")
            }}
            groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            geneVars.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}}}
          reportln("   sampCounts..." + getDateAndTimeString,"debug");
          
          val sampCount_patho = groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              lvl == "PATHO"
            }}.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}
          val sampCount_lpath = groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              lvl == "PATHO" || lvl == "LPATH";
            }}.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}
          val sampCount_patho_canon = groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              lvl_CANON == "PATHO"
            }}.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}
          val sampCount_lpath_canon = groups.map{ grp => {
            val grpSS = groupToSampleMap(grp);
            variantInfo.filter{ case (geneSet,(lvl,lvl_CANON),ss) => {
              lvl_CANON == "PATHO" || lvl_CANON == "LPATH";
            }}.flatMap{ case (geneSet,(lvl,lvl_CANON),ss) => {
              ss
            }}.toSet.filter(s => { grpSS.contains(s) }).size;
          }}
          reportln("   done with all samp counts..." + getDateAndTimeString,"debug");
          
          out.write("geneID");
          groups.map{ grp => {
            out.write("\t"+Vector("numSamp PATHO","numSamp PATHO CANON","numSamp LPATH","numSamp LPATH CANON").map(grp + _).mkString("\t"));
          }}
          out.write("\n");

          allGenes.zipWithIndex.foreach{ case (gene,i) => {
            out.write(gene);
            groups.zipWithIndex.foreach{ case (grp,grpIdx) => {
             out.write("\t"+
                 geneSampCount_patho(i)(grpIdx)+"\t"+
                 geneSampCount_patho_canon(i)(grpIdx)+"\t"+
                 geneSampCount_lpath(i)(grpIdx)+"\t"+
                 geneSampCount_lpath_canon(i)(grpIdx));
            }}
            out.write("\n");
          }}
          
          out.write("ANY_GENE");
          
          groups.zipWithIndex.foreach{ case (grp,grpIdx) => {
             out.write("\t"+
                 sampCount_patho(grpIdx)+"\t"+
                 sampCount_patho_canon(grpIdx)+"\t"+
                 sampCount_lpath(grpIdx)+"\t"+
                 sampCount_lpath_canon(grpIdx));
          }}
          out.write("\n");
          
          out.close();
          reportln("   done with SummaryACMGWalker. " + getDateAndTimeString,"debug");
        }),newHeader);
      }
    
  }
  
  def getAttributeSmart(vc : VariantContext, code : String) : Seq[String] = {
    val attr = vc.getAttributeAsList(code);
    if(attr.length == 1 && attr.head == "."){
      return Seq[String]();
    } else if(attr.length == 0){
      return Seq[String]();
    } else {
      return attr.map(a => a.toString()).toSeq;
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
                 rmskFile : Option[String],
                 refSeqFile : Option[String],
                 inSilicoParams : Option[String] = None,
                 inSilicoMinCt : Option[Int] = None//,
                 //domainSummaryFile : Option[String] = None
                ) extends internalUtils.VcfTool.VCFWalker {
    
    reportln("Creating AssessACMGWalker() ["+stdUtils.getDateAndTimeString+"]","note")
    
    val inSilicoKeysOpt : Option[Seq[(String,Set[String],Set[String])]] = inSilicoParams match {
      case Some(isk) => {
        Some(isk.split(",").toSeq.map(s => {
          val cells = s.split(":");
          if(cells.length != 3){
            error("Malformed input parameter: each comma-delimited element in inSilicoKeys must have 3 colon-delimited parts!");
          }
          (cells(0),cells(1).split("\\|").toSet,cells(2).split("\\|").toSet);
        }))
      }
      case None => None;
    }

    val inSilicoMin = if(inSilicoMinCt.isEmpty){
      if(inSilicoKeysOpt.isEmpty) 0;
      else inSilicoKeysOpt.get.size();
    } else inSilicoMinCt.get;
    
    reportln("Reading txToGene map... ["+stdUtils.getDateAndTimeString+"]","debug");
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
    reportln("Reading geneToTX map... ["+stdUtils.getDateAndTimeString+"]","debug");

    val geneToTx : (String => Set[String]) = txToGeneFile match {
      case Some(f) => {
        val geneToTxMap = scala.collection.mutable.AnyRefMap[String,Set[String]]().withDefault(x => Set[String]());
        getLinesSmartUnzip(f).foreach(line => {
          val cells = line.split("\t");
          geneToTxMap(cells(1)) = geneToTxMap(cells(1)) + cells(0);
        })
        ((s : String) => geneToTxMap.getOrElse(s,Set[String]()));
      }
      case None => {
        ((s : String) => Set[String](s));
      }
    }
    
    reportln("Reading refSeq map... ["+stdUtils.getDateAndTimeString+"]","debug");
    
    val isRefSeq : (String => Boolean) = refSeqFile match {
      case Some(f) => {
        val lines = getLinesSmartUnzip(f);
        val table = getTableFromLines(lines,colNames = Seq("transcript"), errorName = "File "+f);
        var refSeqSet : Set[String] = table.map(tableCells => {
          val tx = tableCells(0);
          tx
        }).toSet
        reportln("   found: "+refSeqSet.size+" RefSeq transcripts.","debug");
        ((s : String) => {
          refSeqSet.contains(s);
        })
      }
      case None => {
        ((s : String) => {
          false;
        })
      }
    }
    /*
    val isRefSeq : (String => Boolean) = refSeqFile match {
      case Some(f) => {
        val lines = getLinesSmartUnzip(f);
        val rawHeaderLine = lines.next();
        val headerCells = if(rawHeaderLine.charAt(0) == '#'){
          rawHeaderLine.tail.split("\t");
        } else {
          rawHeaderLine.split("\t");
        }
        if((! headerCells.contains("name")) || (! headerCells.contains("isRefSeq"))){
          error("FATAL INPUT ERROR: refSeqFile must have at least 2 labelled columns: name and isRefSeq.");
        }
        val nameCol = headerCells.indexOf("name");
        val irsCol = headerCells.indexOf("isRefSeq");
        
        var refSeqSet = Set[String]();
        
        while(lines.hasNext){
          val cells = lines.next().split("\t");
          if(cells(irsCol) == "1"){
            refSeqSet = refSeqSet + cells(nameCol);
          }
        }
        
        reportln("   found: "+refSeqSet.size+" RefSeq transcripts.","debug");
        
        ((s : String) => {
          refSeqSet.contains(s);
        })
      }
      case None => {
        ((s : String) => {
          false;
        })
      }
    }*/
    
    
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
        val tolCol = headerCells.indexWhere(headerCell => headerCell == "MIStolerant");
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
    
    val geneIsMisInsensitive : (String => Boolean) = ((s : String) => { ! geneIsMisSensitive(s); })
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
        val domainIdCol = if(headerCells.contains("domainUID")) headerCells.indexWhere(headerCell => headerCell == "domainUID");
                          else                                  headerCells.indexWhere(headerCell => headerCell == "domainID");
        
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
    
    /*def calcDomainSummary(domainSummaryFile : String, 
                          locusArray : GenomicArrayOfSets[String],
                          locusDomains : ((commonSeqUtils.GenomicInterval) => Set[String]),
                          clinVarVariantSet : scala.collection.Map[String,(String,Int,String)],
                          clinVarVariants : scala.collection.Map[String,Set[internalUtils.TXUtil.pVariantInfo]]){
      
    }*/
    
    def walkVCF(vcIter : Iterator[VariantContext], vcfHeader : VCFHeader, verbose : Boolean = true) : (Iterator[VariantContext],VCFHeader) = {
      //val (clinVarVariantSet,clinVarVariants) : (scala.collection.Map[String,String],scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]]) = getClinVarVariants(infile=clinVarVcf,chromList =chromList, vcfCodes = vcfCodes);
      val (clinVarVariantSet,clinVarVariants) : (scala.collection.Map[String,(String,Int,String)],scala.collection.Map[String,Set[internalUtils.TXUtil.pVariantInfo]]) = getFullClinVarVariants(infile=clinVarVcf,chromList =chromList, vcfCodes = vcfCodes);
      
      //if(domainSummaryFile.isDefined){
        
      //}
      
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

            new VCFInfoHeaderLine(vcfCodes.assess_geneList, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "A list of genes that this variant is on or near."),
            new VCFInfoHeaderLine(vcfCodes.assess_geneTxList, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "For each gene that this variant is on or near, a list of all transcripts belonging to that gene (this is all TX, including TX that this variant is not on or near)."),
            new VCFInfoHeaderLine(vcfCodes.assess_geneLofTxRatio, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "For each gene that this variant overlaps with, the number of transcripts that lose function due to this variant slash the total number of transcripts belonging to that gene."),
            new VCFInfoHeaderLine(vcfCodes.assess_geneMisTxRatio, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "For each gene that this variant overlaps with, the number of transcripts that contain coding missense changes due to this variant slash the total number of transcripts belonging to that gene."),
            
            new VCFInfoHeaderLine(vcfCodes.assess_refSeqLOF, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "List of genes in which the Canonical isoform is loss-of-function due to this variant."),
            new VCFInfoHeaderLine(vcfCodes.assess_refSeqMis, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "List of genes in which the Canonical isoform contains a missense due to this variant"),
            new VCFInfoHeaderLine(vcfCodes.assess_refSeqKnown, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "For each gene that this variant overlaps with, how many canonical transcripts were found (should always be 1 for all genes)."),
            
            new VCFInfoHeaderLine(vcfCodes.assess_LofTX, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "A list of transcripts for which this variant causes loss-of-function (fs, early-stop, start-loss)"),
            new VCFInfoHeaderLine(vcfCodes.assess_MisTX, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "A list of transcripts for which this variant causes a missense."),
            
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
            new VCFInfoHeaderLine(vcfCodes.assess_BP1, 1, VCFHeaderLineType.Integer,      "Missense variant in gene that is NOT missense-sensitive (less than 10pct of pathogenic variants are missense)"),
            new VCFInfoHeaderLine(vcfCodes.assess_BP3, 1, VCFHeaderLineType.Integer,      "In-frame indels in a repetitive region that does NOT overlap with any known domain."),
            new VCFInfoHeaderLine(vcfCodes.assess_BP7, 1, VCFHeaderLineType.Integer,      "Synonymous variant that does NOT intersect with a conserved element region."),
            new VCFInfoHeaderLine(vcfCodes.assess_BA1,    1, VCFHeaderLineType.Integer,   "Allele frequency greater than 5 percent in one or more control dataset ("+ctrlAlleFreqKeys.mkString(",")+")"),
            new VCFInfoHeaderLine(vcfCodes.assess_RATING, 1, VCFHeaderLineType.String,    "ACMG Pathogenicity rating: PATHO - pathogenic. LPATH - likely pathogenic, VUS - variant, unknown significance, LB - likely benign, B - benign."),
            new VCFInfoHeaderLine(vcfCodes.assess_WARNFLAG, 1, VCFHeaderLineType.Integer, "Whether or not there is anything odd about this variant that may require manual inspection."),
            new VCFInfoHeaderLine(vcfCodes.assess_WARNINGS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "List of warnings concerning this variant."),
            
            new VCFInfoHeaderLine(vcfCodes.assess_PVS1_CANON, 1, VCFHeaderLineType.Integer,    "(Considering the canonical TX only) Loss-of-function variant in gene that is LOF-sensitive. Loss-of-function is defined as one of:"+
                                                                                         "stop-gain, frameshift, total-gene-indel, or splice junction indel. "+
                                                                                         "A gene is defined as LOF-sensitive if at least 10pct of pathogenic "+
                                                                                         "clinvar variants are LOF-type variants."),
            new VCFInfoHeaderLine(vcfCodes.assess_PS1_CANON, 1, VCFHeaderLineType.Integer,"(Considering the canonical TX only) Variant has the same amino acid change as a pathogenic variant from ClinVar."),
            new VCFInfoHeaderLine(vcfCodes.assess_PP2_CANON, 1, VCFHeaderLineType.Integer,"(Considering the canonical TX only) Missense variant in gene that is missense-sensitive (missense variants are at least 10pct of known pathogenic variants in clinvar)."),
            new VCFInfoHeaderLine(vcfCodes.assess_BP1_CANON, 1, VCFHeaderLineType.Integer,"(Considering the canonical TX only) Missense variant in gene that is NOT missense-sensitive (less than 10pct of pathogenic variants are missense)"),
            new VCFInfoHeaderLine(vcfCodes.assess_BP3_CANON, 1, VCFHeaderLineType.Integer,"(Considering the canonical TX only) In-frame indels in a repetitive region that does NOT overlap with any known domain."),
            new VCFInfoHeaderLine(vcfCodes.assess_PM4_CANON, 1, VCFHeaderLineType.Integer,      "(Considering the canonical TX only) Protein length change in nonrepeat region OR stop-loss variant."),
            new VCFInfoHeaderLine(vcfCodes.assess_PM5_CANON, 1, VCFHeaderLineType.Integer,      "(Considering the canonical TX only) Novel missense variant change at amino acid where a different amino acid change is known pathogenic in ClinVar."),
            new VCFInfoHeaderLine(vcfCodes.assess_BP7_CANON, 1, VCFHeaderLineType.Integer,      "(Considering the canonical TX only) Synonymous variant that does NOT intersect with a conserved element region."),
            new VCFInfoHeaderLine(vcfCodes.assess_RATING_CANON, 1, VCFHeaderLineType.String,    "(Considering the canonical TX only) ACMG Pathogenicity rating: PATHO - pathogenic. LPATH - likely pathogenic, VUS - variant, unknown significance, LB - likely benign, B - benign."),




            //new VCFInfoHeaderLine(vcfCodes.geneIDs, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,    "List of gene symbols that this variant overlaps with."),
        
            new VCFInfoHeaderLine(vcfCodes.assess_exactMatchInfo, 1, VCFHeaderLineType.String,   "RS number and clinical significance level referring to any ClinVar variant that matches this variant exactly (or blank if there is no matching variant). Clinical significance levels are summarized across all reports. If a variant is reported as both likely benign and benign, it will be marked benign, and likewise with pathogenic. If a variant is marked both pathogenic or likely pathogenic and benign or likely benign, it is marked uncertain. 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other"),
            new VCFInfoHeaderLine(vcfCodes.assess_aminoMatchInfo, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "List of RS numbers and clinical significance level referring to any ClinVar variant(s) that match the amino acid change of this variant  (or blank if there is no matching variant). See "+vcfCodes.assess_exactMatchInfo+" for more info on how the clnsig rating is calculated."),
            new VCFInfoHeaderLine(vcfCodes.assess_nearMatchInfo,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "List of RS numbers and clinical significance level referring to any ClinVar variant(s) that alter the same amino acid position, but change to a different amino acid  (or blank if there is no matching variant). See "+vcfCodes.assess_exactMatchInfo+" for more info on how the clnsig rating is calculated."),

            new VCFInfoHeaderLine(vcfCodes.assess_pathoExactMatchRS, 1, VCFHeaderLineType.String,   "RS number referring to pathogenic ClinVar variant that matches this variant exactly (or blank if there is no matching variant)."),
            new VCFInfoHeaderLine(vcfCodes.assess_pathoAminoMatchRS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "List of RS numbers referring to pathogenic ClinVar variant(s) that match the amino acid change of this variant  (or blank if there is no matching variant)."),
            new VCFInfoHeaderLine(vcfCodes.assess_pathoNearMatchRS,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "List of RS numbers referring to pathogenic ClinVar variant(s) that alter the same amino acid position, but change to a different amino acid  (or blank if there is no matching variant)."),
            
            new VCFInfoHeaderLine(vcfCodes.assess_pathoAminoMatchRS_CANON, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "(Considering the canonical TX only) List of RS numbers referring to pathogenic ClinVar variant(s) that match the amino acid change of this variant  (or blank if there is no matching variant)."),
            new VCFInfoHeaderLine(vcfCodes.assess_pathoNearMatchRS_CANON,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "(Considering the canonical TX only) List of RS numbers referring to pathogenic ClinVar variant(s) that alter the same amino acid position, but change to a different amino acid  (or blank if there is no matching variant)."),
            new VCFInfoHeaderLine(vcfCodes.assess_aminoMatchInfo_CANON, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "(Considering the canonical TX only) List of RS numbers and clinical significance level referring to any ClinVar variant(s) that match the amino acid change of this variant  (or blank if there is no matching variant). See "+vcfCodes.assess_exactMatchInfo+" for more info on how the clnsig rating is calculated."),
            new VCFInfoHeaderLine(vcfCodes.assess_nearMatchInfo_CANON,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "(Considering the canonical TX only) List of RS numbers and clinical significance level referring to any ClinVar variant(s) that alter the same amino acid position, but change to a different amino acid  (or blank if there is no matching variant). See "+vcfCodes.assess_exactMatchInfo+" for more info on how the clnsig rating is calculated.")


            
          ) ++ (
            if(inSilicoKeysOpt.isDefined){
              List[VCFInfoHeaderLine](
                new VCFInfoHeaderLine(vcfCodes.assess_inSilicoSummary,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,   "For each of the following in silico algorithms: ("+inSilicoKeysOpt.get.map(_._1).mkString(",")+") the final determination of whether it was considered \"Damaging\", \"Benign\", or neither (\".\")"),
                new VCFInfoHeaderLine(vcfCodes.assess_PP3, 1, VCFHeaderLineType.Integer,      "Predicted to be damaging by at least "+inSilicoMin+" of the following in silico prediction algorithms: "+inSilicoKeysOpt.get.map(_._1).mkString(",")),
                new VCFInfoHeaderLine(vcfCodes.assess_BP4, 1, VCFHeaderLineType.Integer,      "Predicted to be benign by at least "+inSilicoMin+" of the following in silico prediction algorithms: "+inSilicoKeysOpt.get.map(_._1).mkString(","))
              )
            } else {
              List[VCFInfoHeaderLine]()
            }
          )
            
      val newHeader = internalUtils.VcfTool.addHeaderLines(vcfHeader,newHeaderLines);
      reportln("Walking input VCF...","note")
      return ( 
      vcIter.map(v => {
        assessVariant(v=v,
            clinVarVariants=clinVarVariants,
            txToGene = txToGene, geneToTx = geneToTx,
            geneIsLofSensitive = geneIsLofSensitive,
            geneIsMisSensitive = geneIsMisSensitive,
            geneIsMisInsensitive = geneIsMisInsensitive,
            isRefSeq = isRefSeq,
            locusIsRepetitive = locusIsRepetitive,
            locusIsConserved = locusIsConserved,
            locusIsHotspot = locusIsHotspot,
            ctrlAlleFreqKeys = ctrlAlleFreqKeys,
            inSilicoKeysOpt = inSilicoKeysOpt,
            inSilicoMin = inSilicoMin,
            BA1_AF = BA1_AF,
            PM2_AF = PM2_AF,
            clinVarVariantSet = clinVarVariantSet,
            vcfCodes = vcfCodes)
       }), newHeader );
    }
    reportln("AssessACMGWalker() Created... ["+stdUtils.getDateAndTimeString+"]","note")
  }
  
  def assessVariant(v : VariantContext, 
 //                   clinVarVariants : scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]],
                    clinVarVariants : scala.collection.Map[String,Set[internalUtils.TXUtil.pVariantInfo]],
                    txToGene : (String => String) = ((s : String) => s),
                    geneToTx : (String => Set[String]),
                    geneIsLofSensitive : (String => Boolean) = ((s : String) => true),
                    geneIsMisSensitive : (String => Boolean) = ((s : String) => true),
                    geneIsMisInsensitive : (String => Boolean) = ((s : String) => false),
                    isRefSeq : (String => Boolean),
                    locusIsRepetitive : (commonSeqUtils.GenomicInterval => Boolean) ,
                    locusIsConserved  : (commonSeqUtils.GenomicInterval => Boolean) ,
                    locusIsHotspot    : (commonSeqUtils.GenomicInterval => Boolean) ,
                    ctrlAlleFreqKeys : Seq[String] = Seq("1KG_AF","ESP_EA_AF","ExAC_ALL","SWH_AF_GRP_CTRL"),
                    
                    inSilicoKeysOpt : Option[Seq[(String,Set[String],Set[String])]],
                    inSilicoMin : Int, //-1 == ALL
                    //inSilicoToleratedCt : Int = 1,
                    
                    BA1_AF : Double = 0.05,
                    PM2_AF : Double = 0.0001,
                    clinVarVariantSet : scala.collection.Map[String,(String,Int,String)],
                    vcfCodes : VCFAnnoCodes = VCFAnnoCodes()
                    
                    //clinVarVariants?
                    ) : VariantContext = {
    
    var vb = new htsjdk.variant.variantcontext.VariantContextBuilder(v);
    val crits = new ACMGCritSet();
    val canonCrits = new ACMGCritSet();

    var problemList = Set[String]();
    
    //val inSilicoMin = if(inSilicoMinCt == -1) inSilicoKeys.size else inSilicoMinCt;
    
    val chrom = v.getContig();
    val pos = v.getStart();
    
    val refAlle = v.getReference();
    val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a)).zipWithIndex.filter{case (alt,altIdx) => { alt.getBaseString() != "*" }}

    val txList = v.getAttributeAsList(vcfCodes.txList_TAG).toVector.map(_.toString).filter(_ != ".");
    
    if(txList.length > 0){
    val geneList = txList.map(x => txToGene(x));
    //vb = vb.attribute(vcfCodes.geneIDs , geneList.mkString(","));
    
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
      val canonCombo = combo.filter{ case (g,tx,info,c,i) => isRefSeq(tx) }
      
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
      canonCrits.addCrit(new ACMGCrit("PM",1,isHotspot,Seq[String]()))
      
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
                  "   "+"ref["+refAlle.getBaseString()+"],ALTs=["+altAlleles.map(_._1.getBaseString()).mkString(",")+"]"+key+"=["+v.getAttributeAsList(key).map(_.toString()).mkString(",")+"]"+
                  (if(internalUtils.optionHolder.OPTION_DEBUGMODE) "\n   "+v.toStringWithoutGenotypes() else ""),
                  "POPAF_FORMAT_WARNING",25);
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
      canonCrits.addCrit(new ACMGCrit("BA",1,maxAF >= BA1_AF,Seq[String]()));
      canonCrits.addCrit(new ACMGCrit("PM",2,maxAF <= PM2_AF,Seq[String]()));
      
      //***************************genes:
      
      //val geneList = combo.map{case (g,tx,info,c,i) => g};
      val geneSet = geneList.toSet.toVector.sorted
      vb =  vb.attribute(vcfCodes.assess_geneList , geneSet.padTo(1,".").mkString(","));
      vb =  vb.attribute(vcfCodes.assess_geneTxList , geneSet.map((g)  => { geneToTx(g).toVector.sorted.mkString("|") }).padTo(1,".").mkString(","));
      vb =  vb.attribute(vcfCodes.assess_refSeqKnown, geneSet.map( g => g + ":" + geneToTx(g).count(tx => isRefSeq(tx))).padTo(1,".").mkString(",") );
      
      if(geneSet.size > 1){
        problemList = problemList + "MG|MultipleGene";
      }
      
      val numVariantTypes = Seq[Boolean]( 
        (combo.exists{case (g,tx,info,c,i) => {
          info.severityType == "LLOF" || info.subType == "START-LOSS" || info.subType == "splice"
        }}),
        (combo.exists{case (g,tx,info,c,i) => {
          info.subType == "insAA" || info.subType == "delAA" || info.subType == "indelAA"
        }}),
        (combo.exists{case (g,tx,info,c,i) => {
          info.subType == "swapAA"
        }}),
        (combo.exists{case (g,tx,info,c,i) => {
          info.subType == "cds-synon"
        }})
      ).count(b => b);
      
      if(numVariantTypes > 1){
        problemList = problemList + "MV|MultipleVariantTypes";
      }
      
      //***************************LOF:
      val LofTX = combo.filter{case (g,tx,info,c,i) => {
        info.severityType == "LLOF" || info.subType == "START-LOSS" || info.subType == "splice"
      }}
      
      val LofGenes = LofTX.map{case (g,tx,info,c,i) => {g}}.toSet.toList.sorted;

      vb =  vb.attribute(vcfCodes.assess_LofGenes , LofGenes.padTo(1,".").mkString(","));
      val LofTxList = LofTX.map{case (g,tx,info,c,i) => tx};
      vb =  vb.attribute(vcfCodes.assess_LofTX , LofTxList.padTo(1,".").mkString(","))
      
      val geneLofTx = geneSet.map(g => {
        val txSet = geneToTx(g);
        val txLofCt = LofTxList.count{tx => { txSet.contains(tx) }};
        (g,txLofCt, txSet.size)
      });
      vb =  vb.attribute(vcfCodes.assess_geneLofTxRatio , geneLofTx.map{case (g,a,b) => g+":"+a+"/"+b}.padTo(1,".").mkString(","))
      
      val refSeqLof = geneSet.filter{g => {
        val txSet = geneToTx(g);
        val txLof = LofTxList.filter(tx => {txSet.contains(tx)});
        txLof.exists(tx => {
          isRefSeq(tx);
        })
      }}
      vb =  vb.attribute(vcfCodes.assess_refSeqLOF , refSeqLof.padTo(1,".").mkString(","));
      
      val pvs1flag = LofGenes.exists(g => {
        geneIsLofSensitive(g);
      })
      val pvs1flag_canon = refSeqLof.exists{g => {
        geneIsLofSensitive(g);
      }}
      vb =  vb.attribute(vcfCodes.assess_PVS1 , if(pvs1flag) "1" else "0");
      crits.addCrit(new ACMGCrit("PVS",1,pvs1flag,Seq[String]()));
      canonCrits.addCrit(new ACMGCrit("PVS",1,pvs1flag_canon,Seq[String]()));
      
      vb =  vb.attribute(vcfCodes.assess_PVS1_CANON , if(pvs1flag_canon) "1" else "0");
      
      //***************************Mis:
      val MisTX = combo.filter{case (g,tx,info,c,i) => {
        info.subType == "swapAA";
      }}
      val MisGenes = MisTX.map{case (g,tx,info,c,i) => {
        g
      }}.toSet.toList.sorted;
      
      val MisTxList = MisTX.map{case (g,tx,info,c,i) => tx};
      vb =  vb.attribute(vcfCodes.assess_MisGenes , MisGenes.padTo(1,".").mkString(","));
      vb =  vb.attribute(vcfCodes.assess_MisTX , MisTxList.padTo(1,".").mkString(","));
      
      val geneMisTx = geneSet.map(g => {
        val txSet = geneToTx(g);
        val txMisCt = MisTxList.count{tx => { txSet.contains(tx) }};
        (g,txMisCt, txSet.size)
      });
      vb =  vb.attribute(vcfCodes.assess_geneMisTxRatio , geneMisTx.map{case (g,a,b) => g+":"+a+"/"+b}.padTo(1,".").mkString(","))
      
      val refSeqMis = geneSet.filter{g => {
        val txSet = geneToTx(g);
        val txMis = MisTxList.filter(tx => {txSet.contains(tx)});
        txMis.exists(tx => {
          isRefSeq(tx);
        })
      }}
      vb =  vb.attribute(vcfCodes.assess_refSeqMis , refSeqMis.padTo(1,".").mkString(","));
      
      val pp2flag = MisGenes.exists(g => {
        geneIsMisSensitive(g);
      });
      val bp1flag = (! pvs1flag) && (MisGenes.length > 0) && MisGenes.forall(g => {
        geneIsMisInsensitive(g);
      })
      
      val pp2flag_canon = refSeqMis.exists(g => {
        geneIsMisSensitive(g);
      });
      val bp1flag_canon = (! pvs1flag_canon) && (refSeqMis.length > 0) && refSeqMis.forall(g => {
        geneIsMisInsensitive(g);
      })
      
      vb =  vb.attribute(vcfCodes.assess_PP2 , if(pp2flag) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_BP1 , if(bp1flag) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_PP2_CANON , if(pp2flag_canon) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_BP1_CANON , if(bp1flag_canon) "1" else "0");
       
      crits.addCrit(new ACMGCrit("PP",2,pp2flag,Seq[String]()));
      crits.addCrit(new ACMGCrit("BP",1,bp1flag,Seq[String]()));
      canonCrits.addCrit(new ACMGCrit("PP",2,pp2flag_canon,Seq[String]()));
      canonCrits.addCrit(new ACMGCrit("BP",1,bp1flag_canon,Seq[String]()));
      
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

      val pm4flag_canon = {
        ((! isRepetitive) && ({
           canonCombo.exists{case (g,tx,info,c,i) => {
             info.subType == "insAA" || info.subType == "delAA" || info.subType == "indelAA"
           }}
        })) || canonCombo.exists{case (g,tx,info,c,i) => {
             info.pType == "STOP-LOSS"
           }}
      }
      vb =  vb.attribute(vcfCodes.assess_PM4,if(pm4flag) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_PM4_CANON,if(pm4flag_canon) "1" else "0");
      
      val bp3flag = {
        (isRepetitive) && (! isHotspot) && ({
           combo.exists{case (g,tx,info,c,i) => {
             info.subType == "insAA" || info.subType == "delAA" || info.subType == "indelAA"
           }}
        })
      }
      val bp3flag_canon = {
        (isRepetitive) && (! isHotspot) && ({
           canonCombo.exists{case (g,tx,info,c,i) => {
             info.subType == "insAA" || info.subType == "delAA" || info.subType == "indelAA"
           }}
        })
      }
      vb =  vb.attribute(vcfCodes.assess_BP3,if(bp3flag) "1" else "0");
      vb =  vb.attribute(vcfCodes.assess_BP3_CANON,if(bp3flag_canon) "1" else "0");
      
      crits.addCrit(new ACMGCrit("BP",3,bp3flag,Seq[String]()));
      crits.addCrit(new ACMGCrit("PM",4,pm4flag,Seq[String]()));
      canonCrits.addCrit(new ACMGCrit("BP",3,bp3flag_canon,Seq[String]()));
      canonCrits.addCrit(new ACMGCrit("PM",4,pm4flag_canon,Seq[String]()));
      
      //******************************* ClinVar matching:
      
      var ps1 = ACMGCrit("PS",1,false,Seq[String]());
      var pm5 = ACMGCrit("PM",5,false,Seq[String]());
      var ps1RS = Set[String]();
      var pm5RS = Set[String](); 
      val variantString = v.getContig() + ":"+v.getAttributeAsString(vcfCodes.vMutG_TAG,"?!?");
      val (exactAllRS,exactAllSig,exactAllRawSig) = clinVarVariantSet.getOrElse(variantString,("",1,""))
      val exactRS = if(exactAllSig == 4 || exactAllSig == 5){
        exactAllRS
      } else {
        ""
      }
      
      var aminoMatchInfo = Set[(String,Int,String)]();
      var nearMatchInfo = Set[(String,Int,String)]();
      
      var ps1RS_canon = Set[String]();
      var pm5RS_canon = Set[String]();
      var ps1_canon = ACMGCrit("PS",1,false,Seq[String]());
      var pm5_canon = ACMGCrit("PM",5,false,Seq[String]());
      var aminoMatchInfo_canon = Set[(String,Int,String)]();
      var nearMatchInfo_canon = Set[(String,Int,String)]();
      
      combo.filter{case (g,tx,info,c,i) => {
        info.severityType == "NONSYNON"
      }}.foreach{case (g,tx,info,c,i) => {
        val txvar = clinVarVariants(tx) //.filter(cvinfo => cvinfo.isPatho)
        val exactMatchSet = txvar.filter{cvinfo => {cvinfo.start == info.start && cvinfo.subType == info.subType && cvinfo.pvar == info.pvar}};
        val newExactMatchInfo = exactMatchSet.map(cvinfo => ((cvinfo.ID,cvinfo.CLNSIG,cvinfo.RAWCLNSIG))).toSet
        aminoMatchInfo = aminoMatchInfo ++ newExactMatchInfo
        if(isRefSeq(tx)) aminoMatchInfo_canon = aminoMatchInfo_canon ++ newExactMatchInfo;
        val exactPathoMatch = exactMatchSet.filter(cvinfo => cvinfo.isPatho);
        
        if(exactPathoMatch.size > 0){
          ps1 = ps1.merge(ACMGCrit("PS",1,true,Seq[String](tx + ":" + info.pvar)));
          warning("Adding aminoMatch: rsnums:[\""+exactPathoMatch.map(_.ID).toSet.mkString("\",\"")+"\"] \n"+
                  "                          [\""+exactPathoMatch.map{(cvinfo) => { cvinfo.txid+":"+cvinfo.pvar+":"+cvinfo.ID }}.mkString("\",\"")+"\"]","addAminoMatchRS",50);
          ps1RS = ps1RS ++ exactPathoMatch.map(_.ID).toSet
          
          if(isRefSeq(tx)){
            ps1_canon = ps1_canon.merge(ACMGCrit("PS",1,true,Seq[String](tx + ":" + info.pvar)));
            ps1RS_canon = ps1RS_canon ++ exactPathoMatch.map(_.ID).toSet
          }
        }
        if(info.subType == "swapAA"){
          val partialMatch = txvar.filter{(cvinfo) => {cvinfo.subType == "swapAA" && cvinfo.start == info.start && cvinfo.altAA != info.altAA}};
          val partialPathoMatch = partialMatch.filter(cvinfo => cvinfo.isPatho);
          val newNearMatchInfo = partialMatch.map(cvinfo => ((cvinfo.ID,cvinfo.CLNSIG,cvinfo.RAWCLNSIG))).toSet;
          pm5RS = pm5RS ++ partialPathoMatch.map(_.ID).toSet;
          nearMatchInfo = nearMatchInfo ++ newNearMatchInfo;
          if(isRefSeq(tx)){
            pm5RS_canon = pm5RS_canon ++ partialPathoMatch.map(_.ID).toSet;
            nearMatchInfo_canon = nearMatchInfo_canon ++ newNearMatchInfo;
          }
          if(partialPathoMatch.size > 0 && exactPathoMatch.size == 0){
            pm5 = pm5.merge(ACMGCrit("PM",5,true,Seq[String](tx+":"+info.pvar+"Vs"+partialPathoMatch.head.altAA)));
            if(isRefSeq(tx)) pm5_canon = pm5_canon.merge(ACMGCrit("PM",5,true,Seq[String](tx+":"+info.pvar+"Vs"+partialPathoMatch.head.altAA)));
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
      
      if(ps1_canon.flag){
        vb =  vb.attribute(vcfCodes.assess_PS1_CANON,"1");
        vb =  vb.attribute(vcfCodes.assess_PM5_CANON,"0");
      } else if(pm5_canon.flag){
        vb =  vb.attribute(vcfCodes.assess_PS1_CANON,"0");
        vb =  vb.attribute(vcfCodes.assess_PM5_CANON,"1");
      } else {
        vb =  vb.attribute(vcfCodes.assess_PS1_CANON,"0");
        vb =  vb.attribute(vcfCodes.assess_PM5_CANON,"0");
      }
      canonCrits.addCrit(ps1_canon);
      canonCrits.addCrit(pm5_canon);
      
      vb = vb.attribute(vcfCodes.assess_pathoExactMatchRS, if(exactRS == "") "." else "rs"+exactRS);
      vb = vb.attribute(vcfCodes.assess_pathoAminoMatchRS, ps1RS.toVector.sorted.map("rs"+_.toString()).padTo(1,".").mkString(","));
      vb = vb.attribute(vcfCodes.assess_pathoNearMatchRS,  pm5RS.toVector.sorted.map("rs"+_.toString()).padTo(1,".").mkString(","));
      vb = vb.attribute(vcfCodes.assess_pathoAminoMatchRS_CANON, ps1RS_canon.toVector.sorted.map("rs"+_.toString()).padTo(1,".").mkString(","));
      vb = vb.attribute(vcfCodes.assess_pathoNearMatchRS_CANON,  pm5RS_canon.toVector.sorted.map("rs"+_.toString()).padTo(1,".").mkString(","));
      
      //        assess_exactMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_ExactMatch",
      //  assess_aminoMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_AminoMatch",
      //  assess_nearMatchInfo : String = TOP_LEVEL_VCF_TAG+"ACMG_cvInfo_NearMatch",
      
      if(exactAllRS != ""){
        vb = vb.attribute(vcfCodes.assess_exactMatchInfo,  "rs"+exactAllRS+"|"+exactAllSig+"|"+exactAllRawSig);
      } else {
        vb = vb.attribute(vcfCodes.assess_exactMatchInfo,  ".");
      }
      vb = vb.attribute(vcfCodes.assess_aminoMatchInfo, aminoMatchInfo.toVector.sorted.map{case (a,b,c) => "rs"+a + "|" + b+"|"+c}.padTo(1,".").mkString(","));
      vb = vb.attribute(vcfCodes.assess_nearMatchInfo,  nearMatchInfo.toVector.sorted.map{case (a,b,c) => "rs"+a + "|" + b+"|"+c}.padTo(1,".").mkString(","));
      vb = vb.attribute(vcfCodes.assess_aminoMatchInfo_CANON, aminoMatchInfo_canon.toVector.sorted.map{case (a,b,c) => "rs"+a + "|" + b+"|"+c}.padTo(1,".").mkString(","));
      vb = vb.attribute(vcfCodes.assess_nearMatchInfo_CANON,  nearMatchInfo_canon.toVector.sorted.map{case (a,b,c) => "rs"+a + "|" + b+"|"+c}.padTo(1,".").mkString(","));
      
      //******************************* Insilico:
      //inSilicoMin
      //val numSplit = v.getAttributeAsInt(vcfCodes.numSplit_TAG,1);
      val rawAlles = v.getAttributeAsList(vcfCodes.splitAlle_TAG).map(_.toString());
      
      
      /*
      if(inSilicoKeysOpt.isDefined){
        val inSilicoKeys = inSilicoKeysOpt.get;
        
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
                       "Alt alleles (in order) = "+Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a)).map(_.getBaseString()).mkString(",")+", rawAlles="+rawAlles.mkString(",")+", InSilico("+key+")="+xraw.map(_.toString()).mkString(",")+""+
                       (if(internalUtils.optionHolder.OPTION_DEBUGMODE) "\n   "+v.toStringWithoutGenotypes() else ""), "maybeAmbigInSilicoWarning", limit = 100);
              //warning("One known allele, multiple InSilico fields for key: "+key+"\n"+"line:\n   "+v.toString(), "maybeAmbigInSilicoWarning", limit = 25);
              xraw.map(_.toString()).exists(safeSet.contains(_));
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
                    safeSet.contains(xraw(newidx).toString());
                  } else if(xraw.length == 0){
                    false;
                  } else {
                    problemList = problemList + "IS|inSilicoAlleleMismatch2"
                    warning("Multiple known alleles, wrong # of InSilico fields for key: "+key+"\n"+"line:\n   "+
                            (if(internalUtils.optionHolder.OPTION_DEBUGMODE) "\n   "+v.toStringWithoutGenotypes() else ""),
                            "ambigInSilicoWarning", limit = 25);
                    false;
                  }
              }} >= inSilicoMin;
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
        canonCrits.addCrit(new ACMGCrit("PP",3,isBadFlag,Seq[String]()));
        canonCrits.addCrit(new ACMGCrit("BP",4,isOkFlag,Seq[String]()));
      }*/
      
      if(inSilicoKeysOpt.isDefined){
        val inSilicoKeys = inSilicoKeysOpt.get;
        val inSilicoRaw = inSilicoKeys.map{case (key,damSet,safeSet) => {
            v.getAttributeAsList(key).toList;
        }}
        val inSilicoStats = inSilicoRaw.zip(inSilicoKeys).map{ case (isk,(key,damSet,safeSet)) => {
          val fullSet = damSet ++ safeSet;
          val iskfilt = isk.map(k => k.toString()).filter(fullSet.contains(_));
          if(iskfilt.size == 0){
            "."
          } else if(iskfilt.exists(damSet.contains(_))){
            "Damaging"
          } else if(iskfilt.exists(safeSet.contains(_))){
            "Benign"
          } else {
            "."
          }
        }}
        
        val inSilicoBenignCt = inSilicoStats.count(_ == "Benign");
        val inSilicoDamCt = inSilicoStats.count(_ == "Damaging");
        
        val isDamaging = inSilicoDamCt >= inSilicoMin;
        val isBenign   = inSilicoBenignCt >= inSilicoMin;
        
        crits.addCrit(new ACMGCrit("PP",3,isDamaging,Seq[String]()));
        crits.addCrit(new ACMGCrit("BP",4,isBenign,Seq[String]()));
        canonCrits.addCrit(new ACMGCrit("PP",3,isDamaging,Seq[String]()));
        canonCrits.addCrit(new ACMGCrit("BP",4,isBenign,Seq[String]()));
        vb = vb.attribute(vcfCodes.assess_PP3, if(isDamaging) "1" else "0" );
        vb = vb.attribute(vcfCodes.assess_BP4, if(isBenign) "1" else "0" );
        
        vb = vb.attribute(vcfCodes.assess_inSilicoSummary, inSilicoStats.mkString(",") );
      }
      
      
      //******************************* BP7: Synonymous w/ no splice impact, not highly conserved:
      val bp7flag = (! isConserved) && combo.forall{case (g,tx,info,c,i) => {
        info.severityType == "SYNON" || info.severityType == "PSYNON" || info.severityType == "UNK"
      }}
      vb = vb.attribute(vcfCodes.assess_BP7, if(bp7flag) "1" else "0" );
      crits.addCrit(new ACMGCrit("BP",7,bp7flag,Seq[String]()));
      
      val bp7flag_CANON = (! isConserved) && canonCombo.forall{case (g,tx,info,c,i) => {
        info.severityType == "SYNON" || info.severityType == "PSYNON" || info.severityType == "UNK"
      }}
      vb = vb.attribute(vcfCodes.assess_BP7_CANON, if(bp7flag_CANON) "1" else "0" );
      canonCrits.addCrit(new ACMGCrit("BP",7,bp7flag_CANON,Seq[String]()));
      
      //*******************************
      
      for((critLvl, critNum, critTAG) <- VcfTool.UNAUTOMATED_ACMG_PARAMS){
        if(v.hasAttribute(critTAG)){
          crits.addCrit(new ACMGCrit(critLvl,critNum,v.getAttributeAsString(critTAG,"") == "1",Seq[String]()));
          canonCrits.addCrit(new ACMGCrit(critLvl,critNum,v.getAttributeAsString(critTAG,"") == "1",Seq[String]()));
        }
      }
      
      val rating = CalcACMGVar.getACMGPathogenicityRating(crits.getCritSet);
      val rating_CANON = CalcACMGVar.getACMGPathogenicityRating(canonCrits.getCritSet);
      
      vb = vb.attribute(vcfCodes.assess_RATING, rating );
      vb = vb.attribute(vcfCodes.assess_RATING_CANON, rating_CANON );

    //vb =  vb.attribute(vcfCodes.assess_PP3,if(inSilicoFlag,"1","0"));
      //val inSilico
      
      
      //BP3:
      
      //crits;
    //}}.toVector;
    
    //return (out,vb.make());
    }
    
    vb = vb.attribute(vcfCodes.assess_WARNFLAG, if(problemList.isEmpty) "0" else "1" );
    vb = vb.attribute(vcfCodes.assess_WARNINGS, problemList.toVector.sorted.padTo(1,".").mkString(",") );
    
    return vb.make();
  }
  
  def getClinVarVariants(infile : String, chromList : Option[List[String]], vcfCodes : VCFAnnoCodes = VCFAnnoCodes()) : (scala.collection.Map[String,String],scala.collection.Map[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]]) = {
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(infile, 
                                                                  chromList = chromList,
                                                                  vcfCodes = vcfCodes);
    
    val out = new scala.collection.mutable.AnyRefMap[String,Set[(String,internalUtils.TXUtil.pVariantInfo,String)]](((k : String) => Set[(String,internalUtils.TXUtil.pVariantInfo,String)]()));
    var gset = new scala.collection.mutable.AnyRefMap[String,String]();
    
    reportln("Starting VCF read/write...","progress");
    for(v <- vcIter){
      try {
      val rsnum = v.getAttributeAsString("RS","unknownRSNUM") //v.getID();
      val refAlle = v.getReference();
      val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a));
      //val vTypesList = v.getAttributeAsList(vcfCodes.vType_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
      val txList = v.getAttributeAsList(vcfCodes.txList_TAG).toVector.map(_.toString).filter(_ != ".");
      //val vMutPList = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
      val vMutCListRaw = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.filter(_ != ".");
      val vMutCList = vMutCListRaw.map(_.toString.split("\\|").toVector);
      val chrom = v.getContig();
      
      val vMutGList = v.getAttributeAsList(vcfCodes.vMutG_TAG).toVector.filter(_ != ".");

      
      val vMutInfoList = if(txList.length > 0) {
          v.getAttributeAsList(vcfCodes.vMutINFO_TAG).toVector.map( (attrObj) => {
          val attrString = attrObj.toString();
          val attrSplit = attrString.split("\\|").toVector
          //reportln("vMutInfoListDebug: attrString = \""+attrString+"\"","debug");
          //reportln("vMutInfoListDebug: attrSplit["+attrSplit.length+"] = [\""+attrSplit.mkString("\", \"")+"\"]","debug");
          
          attrSplit.map(x => {
            internalUtils.TXUtil.getPvarInfoFromString(x, ID = rsnum);
          })
        })
      } else {
        Vector();
      }
      val vClnSig = v.getAttributeAsList("CLNSIG").toVector.map(_.toString());
      
      
      if(txList.length > 0) { 
        for((alle,altIdx) <- altAlleles.zipWithIndex){
          
      //vMutGList.foreach(g => {
      //  gset(chrom + ":" + g) = rsnum;
      //})
           
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
               out(tx) = out(tx) + ((vMutC(i),vInfo(i),rsnum));
             }
             gset(chrom + ":" + vMutGList(altIdx)) = rsnum;
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
    
    out.keySet.toVector.slice(0,10).foreach(tx => {
      val varSeq = out(tx);
      reportln("Example ClinVar TX: "+tx,"debug");
      varSeq.slice(0,10).foreach{ case (cvc,info,rsnum) => {
        reportln("      " + info.txid + ":" + info.pvar + ":"+rsnum,"debug");
      }}
    })
    
    return (gset,out);
  }
  
  //output: tx => Set[(HGVDc,pVariantInfo])]
  def getFullClinVarVariants(infile : String, chromList : Option[List[String]], vcfCodes : VCFAnnoCodes = VCFAnnoCodes()) : (scala.collection.Map[String,(String,Int,String)],scala.collection.Map[String,Set[internalUtils.TXUtil.pVariantInfo]]) = {
    val (vcIter,vcfHeader) = internalUtils.VcfTool.getVcfIterator(infile, 
                                                                  chromList = chromList,
                                                                  vcfCodes = vcfCodes);
    
    val out = new scala.collection.mutable.AnyRefMap[String,Set[internalUtils.TXUtil.pVariantInfo]](((k : String) => Set[internalUtils.TXUtil.pVariantInfo]()));
    var gset = new scala.collection.mutable.AnyRefMap[String,(String,Int,String)]();
    
    reportln("Starting VCF read/write...","progress");
    for(v <- vcIter){
      try {
      val rsnum = v.getAttributeAsString("RS","unknownRSNUM") //v.getID();
      val refAlle = v.getReference();
      val altAlleles = Range(0,v.getNAlleles()-1).map((a) => v.getAlternateAllele(a));
      //val vTypesList = v.getAttributeAsList(vcfCodes.vType_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
      val txList = v.getAttributeAsList(vcfCodes.txList_TAG).toVector.map(_.toString).filter(_ != ".");
      //val vMutPList = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.map(_.toString.split(vcfCodes.delims(1)).toVector);
      val vMutCListRaw = v.getAttributeAsList(vcfCodes.vMutP_TAG).toVector.filter(_ != ".");
      val vMutCList = vMutCListRaw.map(_.toString.split("\\|").toVector);
      val chrom = v.getContig();
      
      val vMutGList = v.getAttributeAsList(vcfCodes.vMutG_TAG).toVector.filter(_ != ".");

      val vClnSig = v.getAttributeAsList("CLNSIG").toVector.map(_.toString());
      val clnsig = vClnSig.zipWithIndex.map{case (sigstr,altidx) => {
        val clnSig = sigstr.split("\\|");
        (getSummaryClinSig(clnSig),clnSig.mkString(":"));
      }}

      val vMutInfoList = if(txList.length > 0) {
        v.getAttributeAsList(vcfCodes.vMutINFO_TAG).toVector.zipWithIndex.map{ case (attrObj, altIdx) => {
          val attrString = attrObj.toString();
          val attrSplit = attrString.split("\\|").toVector
          //reportln("vMutInfoListDebug: attrString = \""+attrString+"\"","debug");
          //reportln("vMutInfoListDebug: attrSplit["+attrSplit.length+"] = [\""+attrSplit.mkString("\", \"")+"\"]","debug");
          
          attrSplit.map{ case (x) => {
            //val clnSig = vClnSig(altIdx).split("\\|");
            //0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other                        
            val (cs,csRaw) = clnsig(altIdx) //getSummaryClinSig(clnSig);
            internalUtils.TXUtil.getPvarInfoFromString(x, ID = rsnum,CLNSIG=cs,RAWCLNSIG=csRaw);
          }}
        }}
      } else {
        Vector();
      }
      
      if(txList.length > 0) { 
        for((alle,altIdx) <- altAlleles.zipWithIndex.filter{case (a,i) => { a.getBaseString() != "*" }}){
           //vMutGList.foreach(g => {
           //  gset(chrom + ":" + g) = rsnum;
           //})
           //val vTypes = vTypesList(altIdx);
           //val vMutP = vMutPList(altIdx);
           val vMutC = vMutCList(altIdx);
           val vInfo = vMutInfoList(altIdx);
           
           if(vMutC.length != txList.length){
             reportln("vMutC.length = "+vMutC.length+", txList = "+txList+"\nvMutC = [\""+vMutC.mkString("\",\"")+"\"]","debug");
           }
           //if(hasPatho && (! hasBenign)){
             for(i <- Range(0,txList.length)){
               val tx = txList(i);
               out(tx) = out(tx) + (vInfo(i));
             }
             gset(chrom + ":" + vMutGList(altIdx)) = (rsnum,clnsig(altIdx)._1,clnsig(altIdx)._2);
           //}
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
    
    out.keySet.toVector.slice(0,10).foreach(tx => {
      val varSeq = out(tx);
      reportln("Example ClinVar TX: "+tx,"debug");
      varSeq.slice(0,10).foreach{ (info) => {
        reportln("      " + info.txid + "|" + info.pvar + "|"+info.ID+"|"+info.CLNSIG+"|"+info.RAWCLNSIG,"debug");
      }}
    })
    
    return (gset,out);
  }
  
  def getSummaryClinSig(clnSig : Seq[String]) : Int = {
      var cs = 1;
      
      if( clnSig.contains("1")){
        cs = 1;
      }
      if( clnSig.contains("0") ){
        cs = 0;
      }
      if( clnSig.contains("255")){
       cs = 255;
      }
      if( clnSig.contains("6")){
        cs = 6;
      }
      if( clnSig.contains("7")){
        cs = 7;
      }
      if(clnSig.contains("3")){
        cs = 3;
      }
      if(clnSig.contains("2")){
        cs = 2;
      }
      if(clnSig.contains("4") && (cs == 2 || cs == 3)){
        return 0;
      }
      if(clnSig.contains("4")){
        cs = 4;
      }
      if(clnSig.contains("5") && (cs == 2 || cs == 3)){
        return 0;
      }
      if(clnSig.contains("5")){
        cs = 5;
      }
      
      return cs;
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









