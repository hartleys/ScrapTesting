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
import scala.collection.JavaConverters._
import internalUtils.optionHolder._;
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import internalUtils.commandLineUI._;
import internalUtils.fileUtils._;

import htsjdk.variant._;
import htsjdk.variant.variantcontext._;
import htsjdk.variant.vcf._;

//import javaInternalUtils.TXUtilEnum._;
 

import internalUtils.genomicUtils._;
import internalUtils.commonSeqUtils._;

  /** An annotation container and utility class. Instances of this class contains various information 
   *  regarding one specific transcript and provide a variety of helper functions to 
   *  query this information.
   * 
   *  In general the constructor should not be used externally. It should be created using the 
   *  buildTXUtilsFromAnnotation() function.
   *  
   *  All coordinate positions are 0-based.
   *  
   *  @Param geneID           The geneID for the transcript.
   *  @Param txID             The transcript ID
   *  @Param chrom            The chromosome on which the TX is found.
   *  @Param strand           The strand on which the TX is found. Either '+' or '-'.
   *  @Param gSpans           The genomic coordinates for each exome, in order of ascending genomic-coordinates.
   *  @Param gCdsStart        The earliest genomic coordinate of the CDS (coding) region, including the stop codon.
   *  @Param gCdsEnd          The latest genomic coordinate of the CDS (coding) region, including the stop codon.
   *  @Param seqGS            The sequence of base-pairs for the mRNA, on the genomic strand.
   *  @Param leadingExonSeq   A list of short strings consisting of the genomic bases immediately before each exon.
   *  @Param trailingExonSeq  A list of short strings consisting of the genomic bases immediately after each exon.
   *  @Param gCdsNS           The start and end genomic coords for the CDS region, not counting the stop codon.
   */   
case class TXUtil(geneID : String, txID : String, chrom : String, strand : Char, 
                  gSpans : Vector[(Int,Int)], 
                  gCdsStart : Int, gCdsEnd : Int,
                  seqGS : Vector[Char],
                  leadingExonSeq : Vector[String],
                  trailingExonSeq : Vector[String],
                  gCdsNS : (Int,Int)
                  ){
  
  lazy val gStart = gSpans.head._1;
  lazy val gEnd = gSpans.last._2;
  
  lazy val gExonIDs = if(strand == '+') Range(0,numExons) else Range(0,numExons).reverse;
  
  lazy val gIntronSpans = gSpans.init.zip(gSpans.tail).map{ case ((a,b),(c,d)) => { (b,c) }};
  
  /** The coordinates for each exome on the mRNA, read in the genomic direction.
   *  
   */
  lazy val rSpansGS : Vector[(Int,Int)] = {
      gSpans.tail.foldLeft(Vector[(Int,Int)]((0,gSpans.head._2 - gSpans.head._1))){ case (soFar,curr) => {
        soFar :+ ((soFar.last._2, soFar.last._2 + (curr._2 - curr._1) ));
      }}
    }
  /** The mRNA sequence, on the mRNA strand.
   * 
   */
  lazy val seq : Vector[Char] = if(strand == '+') seqGS else internalUtils.commonSeqUtils.complementSequenceString(seqGS.reverse).toVector;
  /** A short string consisting of the genomic bases immediately before each exon, read on the mRNA strand.
   *  (Note that if the TX is on the reverse strand, the exons will be in reverse order)
   */
  lazy val rLeadingExonSeq : Vector[String] = if(strand == '+') leadingExonSeq else {
      trailingExonSeq.reverse.map(s => {
        internalUtils.commonSeqUtils.complementSequenceString(s.reverse).mkString("");
      })
    }
  /** A short string consisting of the genomic bases immediately following each exon, read on the mRNA strand.
   *  (Note that if the TX is on the reverse strand, the exons will be in reverse order)
   */
  lazy val rTrailingExonSeq : Vector[String] = if(strand == '+') trailingExonSeq else {
      leadingExonSeq.reverse.map(s => {
        internalUtils.commonSeqUtils.complementSequenceString(s.reverse).mkString("");
      })
    }
  lazy val cTrailingSeq : String = seq.drop(rCdsEnd).mkString("") + rTrailingExonSeq.last;
  //lazy val stopCodon : String = cTrailingSeq.substring(0,3);
  
  /** The number of exons in the transcript.
   *  
   */
  lazy val numExons : Int = gSpans.length;
  /** The mRNA coordinates for each exon, read in the mRNA strand direction.
   *  (Note that if the TX is on the reverse strand, the exons will be in reverse order)
   */
  lazy val rSpans : Vector[(Int,Int)] = if(strand == '+'){
      rSpansGS;
    } else {
      rSpansGS.reverse.tail.foldLeft(Vector((0,rSpansGS.last._2 - rSpansGS.last._1))){case (soFar,curr) => {
        soFar :+ ((soFar.last._2, soFar.last._2 + (curr._2 - curr._1) ));
      }}
    }
  lazy val cSpans : Vector[(Int,Int)] = {
      val cdsExons = rSpans.filter{case (s,e) => s < rCdsEnd && rCdsStart < e}
      if(cdsExons.length == 1){
        Vector((rCdsStart,rCdsEnd));
      } else if(cdsExons.length == 2){
        Vector((rCdsStart,cdsExons.head._2), (cdsExons.last._1,rCdsEnd))
      } else {
        Vector((rCdsStart,cdsExons.head._2)) ++ cdsExons.init.tail ++ Vector((cdsExons.last._1,rCdsEnd))
      }
    }
  
  lazy val gSpansCDS : Vector[(Int,Int)] = {
      val cdsExons = gSpans.filter{case (s,e) => s < gCdsEnd && gCdsStart < e}
      if(cdsExons.length == 1){
        Vector((gCdsStart,gCdsEnd));
      } else if(cdsExons.length == 2){
        Vector((gCdsStart,cdsExons.head._2), (cdsExons.last._1,gCdsEnd))
      } else {
        Vector((gCdsStart,cdsExons.head._2)) ++ cdsExons.init.tail ++ Vector((cdsExons.last._1,gCdsEnd))
      }
    }
  //lazy val gSpansCDSRS : Vector[(Int,Int)] = if(strand == '+') gSpansCDS else gSpansCDS.reverse;
  /** A Vector of each exons genomic sequence, on the genomic strand.
   * 
   */
  lazy val exonSeqGS : Vector[String] = {
      rSpansGS.map{case (start,end) =>{
        seqGS.slice(start,end).mkString("")
      }}
    }
  lazy val exonSeq : Vector[String] = {
    rSpans.map{case (start,end) => {
      seq.slice(start,end).mkString("")
    }}
  }
  
  /** The number of genomic bases spanned by the TX.
   * 
   */
  lazy val gLen = gSpans.last._2 - gSpans.head._1;
  /** The length of the mRNA, in base-pairs.
   * 
   */
  lazy val rLen = rSpans.last._2;
  
  /** Converts a genomic coordinate to a coordinate on the mRNA.
   *  
   */
  def convertGtoR(g : Int, strict : Boolean = false) : Option[Int] = {
      val idx = gSpans.indexWhere{ case (gs,ge) => {
        g >= gs && g < ge
      }}
      if(idx == -1) {
        if(strict) reportln("Attempted convertGtoC("+g+") for tx:\n"+toString(),"debug");
        None;
      } else if(strand == '+'){
        val distFromExonStart = g - gSpans(idx)._1;
        Some(rSpansGS(idx)._1 + distFromExonStart);
      } else {
        val distFromExonStart = gSpans(idx)._2 - 1 - g;
        val rIdx = rSpans.length - idx - 1;
        Some(rSpans(rIdx)._1 + distFromExonStart);
      }
    }
  
  /** The mRNA coordinate of the start of the CDS region.
   *  
   */
  lazy val rCdsStart = if(strand == '+') convertGtoR(gCdsStart).get     else convertGtoR(gCdsEnd - 1).get
  /** The mRNA coordinate of the end of the CDS region. Coordinates are 0-based, and thus this is actually the coordinate 1 base after the last coding base.
   * 
   */
  lazy val rCdsEnd =   if(strand == '-') convertGtoR(gCdsStart).get + 1 else convertGtoR(gCdsEnd - 1).get + 1
  
  lazy val rUtrLen5 : Int = rCdsStart;
  lazy val rUtrLen3 : Int = rLen - rCdsEnd;
  
  private lazy val gSpansZipcSpans : Vector[((Int,Int),(Int,Int))] = if(strand == '+') gSpansCDS.zip(cSpans) else gSpansCDS.zip(cSpans.reverse);
  private lazy val cSpansZipgSpans : Vector[((Int,Int),(Int,Int))] = if(strand == '+') cSpans.zip(gSpansCDS) else cSpans.zip(gSpansCDS.reverse);
  /** The coding sequence.
   *  
   */
  lazy val cSeq : Vector[Char] = {
      seq.slice(rCdsStart,rCdsEnd);
    }
  /** The length of the coding sequence.
   * 
   */
  lazy val cLen = rCdsEnd - rCdsStart;
  
  /** A vector of the coding codons.
   * 
   */
  lazy val codons : Vector[String] = internalUtils.commonSeqUtils.getCodonsFromSeq(cSeq)// :+ stopCodon;
  /** A vector of the amino acids, using the standard 1-character abbreviations.
   * 
   */
  lazy val aa : Seq[String] = codons.map(c => internalUtils.commonSeqUtils.getAminoAcidFromCodon(c,0,true));
  lazy val aminos : Seq[String] = codons.view.map(c => internalUtils.commonSeqUtils.getAminoAcidFromCodon(c,1,true));
  /** Whether or not the transcript is "valid".
   *  
   */
  def isValidFullLenTX : Boolean = {
      (cLen >= 6) &&
      (cLen % 3 == 0) &&
      (aa.head == "M") &&
      (aa.last == "*") &&
      (aa.count(a => a == "*") == 1)
    }
  
  /** A list of strings of potential issues found with this transcript.
   *  
   */
  def getIssueList : Vector[String] = {
      var out : Vector[String] = Vector[String]();
      if(cLen % 3 != 0) out = out :+ "cdsLenNotDivByThree";
      if(aa.head != "M") out = out :+ "cdsDoesNotStartWithMet";
      if(aa.last != "*"){
        if(internalUtils.commonSeqUtils.getAminoAcidFromCodon(rTrailingExonSeq.last.substring(0,3),0,false) != "*"){
          //stop codon right after CDS region
        } else {
          out = out :+ "noStopCodon";
        }
      } else {
        //contains stop codon
      }
      if(aa.init.contains("*")) out = out :+ "cdsContainsEarlyStopCodon";
      
      return out;
    }
    
  /** String representation of the information included in the TX.
   *  
   */
  override def toString() : String = {
     "Transcript(geneID="+geneID+", txID="+txID+", chrom="+chrom+", strand="+strand+",\n  gCdsStart="+gCdsStart+", gCdsEnd="+gCdsEnd+", ...)\n"+
     "  rCdsStart="+rCdsStart+", rCdsEnd="+rCdsEnd+"\n"+
     "  gLen="+gLen+", rLen="+rLen+", cLen="+cLen+"\n"+
     "  "+this.gSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     "  "+this.rSpansGS.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     "  "+this.rSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     ""+this.seqGS.mkString("")+"\n"+
     ""+this.rSpansGS.map{case (i,j) => repString(" ",j-i-1) + "|"}.mkString("")+"\n"+
     ""+this.seq.mkString("")+"\n"+
     ""+this.rSpans.map{case (i,j) => repString(" ",j-i-1) + "|"}.mkString("")+"\n"+
     ""+this.cSeq.mkString("")+"\n"+
     //"  "+this.codons.mkString(" ")+"\n"+
     " "+this.aa.mkString("  ")+"\n"+
     ""
    } 
  /** A more extensive representation of the information included in the TX.
   * 
   */
  def toStringVerbose() : String = {
     "Transcript(geneID="+geneID+", txID="+txID+", chrom="+chrom+", strand="+strand+",\n  gCdsStart="+gCdsStart+", gCdsEnd="+gCdsEnd+", ("+gCdsNS._1+","+gCdsNS._2+"), ...)\n"+
     "  rCdsStart="+rCdsStart+", rCdsEnd="+rCdsEnd+"\n"+
     "  gLen="+gLen+", rLen="+rLen+", cLen="+cLen+"\n"+
     "  "+this.gSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     "  "+this.rSpansGS.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     "  "+this.rSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     ""+this.seqGS.mkString("")+"\n"+
     ""+gSpans.map{case (s,e) => Range(s,e)}.flatten.map(i => {
         if(i == gCdsNS._1){"<"} else if(i == gCdsNS._2){">"} else if(gSpans.exists(_._1 == i)){"["} else if(gSpans.exists(_._2 - 1 == i)){"]"} else {" "}
       }).mkString("") +"\n"+
     ""+this.seq.mkString("")+"\n"+
     ""+this.rSpans.map{case (i,j) => repString(" ",j-i-1) + "["}.mkString("")+"\n"+     
     ""+this.cSeq.mkString("")+"\n"+
     //"  "+this.codons.mkString(" ")+"\n"+
     " "+this.aa.mkString("  ")+"\n"+
     ""+repString(" ",cSeq.length)+this.cTrailingSeq+"\n"+   
     ""+this.exonSeqGS.zip(leadingExonSeq).zip(trailingExonSeq).zip(gSpans).map{ case ((((exonSeq,leadSeq),tailSeq),(start,end))) => {
        chrom+":"+start+"-"+end+" "+leadSeq +" "+ exonSeq +" "+ tailSeq
     }}.mkString("\n") + "\n"+
     ""
    }
   
  def toStringShort() : String = {
     "Transcript(geneID="+geneID+", txID="+txID+", chrom="+chrom+", strand="+strand+",\n  gCdsStart="+gCdsStart+", gCdsEnd="+gCdsEnd+", ...)\n"+
     "  "+this.gSpans.map{case (i,j) => "("+i + ","+j+")"}.mkString(" ")+"\n"+
     ""+this.seqGS.mkString("")+"\n"+
     ""+gSpans.map{case (s,e) => Range(s,e)}.flatten.map(i => {
         if(i == gCdsStart){"<"} else if(i == gCdsEnd){">"} else if(gSpans.exists(_._1 == i)){"|"} else {" "}
       }).mkString("")+ "\n"       
    }
           
  private def reverseIfReverse[A](sq : Seq[A]) : Seq[A] = {
    if(strand == '+') sq;
    else sq.reverse;
  }
  
  def convertGtoC(gpos : Int) : Option[Int] = {
    convertGtoR(gpos,false) match {
      case Some(r) => {
        if(r >= rCdsStart && r < rCdsEnd){
          Some(r-rCdsStart)
        } else None;
      }
      case None => None;
    }
  }
  def convertCtoG(c : Int) : Option[Int] = {
    val r = c + rCdsStart;
    cSpansZipgSpans.find{case ((rS,rE),(gS,gE)) => {
      r >= rS  && r < rE
    }} match {
      case Some(((rS,rE),(gS,gE))) => {
        if(strand == '+') return Some((r - rS) + gS)
        else return Some(gE - (r - rS) - 1);
      }
      case None => return None;
    }
  }
  
  def convertCSpanToGSpan(start : Int, end : Int) : Vector[(Int,Int)] = {
    val (rStart,rEnd) = (start + rCdsStart, end + rCdsStart);
    
    var ix = cSpansZipgSpans.filter{case ((rS,rE),(gS,gE)) => {
      rStart < rE && rEnd > rS;
    }}
    if(ix.length == 0){
      return Vector[(Int,Int)]();
    }
    
    if(strand == '+'){
      ix = Vector(((rStart,ix.head._1._2),(ix.head._2._1 +  rStart - ix.head._1._1,ix.head._2._2))) ++ ix.tail;
      ix = ix.init ++ Vector(((ix.last._1._1,rEnd  ),(ix.last._2._1, ix.last._2._2 - (ix.last._1._2 - rEnd)  )));
    } else {
      ix = Vector(((rStart,ix.head._1._2),(ix.head._2._1, ix.head._2._2 - (rStart - ix.head._1._1) ))) ++ ix.tail;
      ix = ix.init ++ Vector(((ix.last._1._1,rEnd  ),(ix.last._2._1 + (ix.last._1._2 - rEnd) , ix.last._2._2)));
    }
    return ix.map{ case ((rS,rE),(gS,gE)) => (gS,gE)}
  }
  
  def getCodonPos(gpos: Int) : Option[Int] = {
    convertGtoC(gpos) match {
      case Some(c) => {
        Some(c / 3);
      }
      case None => None;
    }
  }
  
  def getCPos(g : Int) : String = {
    convertGtoC(g) match {
      case Some(c) => {
        return ""+(c+1);
      }
      case None => {
        convertGtoR(g) match {
          case Some(r) => {
            if(r < rCdsStart){
              return "" + "-" + (rCdsStart - r)
            } else if(r >= rCdsEnd){
              return "" + "*" + (r - rCdsEnd + 1)
            } else {
              //IMPOSSIBLE STATE!
              error("getCPos: Impossible state!");
              return "???";
            }
          }
          case None => {
            if(g < gStart){
              val utrlen = if(strand == '+') rUtrLen5 else rUtrLen3;
              val typeChar = if(strand == '+') "-" else "*";
              "" + typeChar + (utrlen + gStart - g)
            } else if(g >= gEnd){
              //3' genomic end, off the end of the annotated UTR
              val utrlen = if(strand == '+') rUtrLen3 else rUtrLen5;
              val typeChar = if(strand == '+') "*" else "-";
              "" + typeChar + (utrlen + g - gEnd + 1)
            } else {
              //NEAR SPLICE
              val (upCoord,dnCoord) = (gSpans.unzip._1,gSpans.unzip._2.map(_-1));
              val distToNearestSJup = upCoord.map(p => math.abs(g - p)).min;
              val distToNearestSJdn = dnCoord.map(p => math.abs(g - p)).min;
              val distToSJ = math.min(distToNearestSJup,distToNearestSJdn);
              val upIsCloser = if(distToNearestSJup == distToNearestSJdn){
                strand == '-'
              } else {
                distToNearestSJup < distToNearestSJdn
              }
              val typeChar = if(strand == '+' && upIsCloser) "-" else if(strand == '+' && ! upIsCloser) "+" else if(strand == '-' && upIsCloser) "+" else "-";
              val pos = if(upIsCloser) upCoord.find(p => math.abs(g-p) == distToSJ).get;
                        else dnCoord.find(p => math.abs(g-p) == distToSJ).get;
              return getCPos(pos) + typeChar + distToSJ;
            }
          }
        }
      }
    }
  }
  
  def getRPos(g : Int) : String = {
    convertGtoR(g) match {
      case Some(r) => {
        return "" + r;
      }
      case None => {
        if(g < gStart){
          val typeChar = if(strand == '+') "-" else "*";
          "" + typeChar + (gStart - g)
        } else if(g >= gEnd){
          val typeChar = if(strand == '+') "*" else "-";
          "" + typeChar + (g - gEnd + 1)
        } else {
          //NEAR SPLICE
          val (upCoord,dnCoord) = (gSpans.unzip._1,gSpans.unzip._2.map(_-1));
          val distToNearestSJup = upCoord.map(p => math.abs(g - p)).min;
          val distToNearestSJdn = dnCoord.map(p => math.abs(g - p)).min;
          val distToSJ = math.min(distToNearestSJup,distToNearestSJdn);
          val upIsCloser = if(distToNearestSJup == distToNearestSJdn){
            strand == '-'
          } else {
            distToNearestSJup < distToNearestSJdn
          }
          val typeChar = if(strand == '+' && upIsCloser) "-" else if(strand == '+' && ! upIsCloser) "+" else if(strand == '-' && upIsCloser) "+" else "-";
          val pos = if(upIsCloser) upCoord.find(p => math.abs(g-p) == distToSJ).get;
                    else dnCoord.find(p => math.abs(g-p) == distToSJ).get;
          return getRPos(pos) + typeChar + distToSJ;
        }
      }
    }
  }
  
  private def getOffsetFlankSeqFunction(gStartPos : Int) : (Int => Option[Char]) = {
    val offsetSeq = if(strand == '+'){
      gSpans.zipWithIndex.find{case ((s,e),exonIdx) => {
        s <= gStartPos && gStartPos < e;
      }} match {
        case Some(((s,e),exonIdx)) => {
          exonSeqGS(exonIdx).drop(gStartPos - s) + rTrailingExonSeq(exonIdx);
        }
        case None => {
          gSpans.zipWithIndex.reverse.find{case ((s,e),exonIdx) => {
            gStartPos >= e && gStartPos < e + rTrailingExonSeq(exonIdx).length;
          }} match {
            case Some(((s,e),exonIdx)) => {
              rTrailingExonSeq(exonIdx).drop(gStartPos - e);
            }
            case None => "";
          }
        }
      }
    } else {
      gSpans.reverse.zipWithIndex.find{case ((s,e),exonIdx) => {
        s <= gStartPos && gStartPos < e;
      }} match {
        case Some(((s,e),exonIdx)) => {
          exonSeq(exonIdx).drop(e - gStartPos - 1) + rTrailingExonSeq(exonIdx);
        }
        case None => {
          gSpans.reverse.zipWithIndex.reverse.find{case ((s,e),exonIdx) => {
            gStartPos < s && gStartPos >= s - rTrailingExonSeq(exonIdx).length;
          }} match {
            case Some(((s,e),exonIdx)) => {
              rTrailingExonSeq(exonIdx).drop(s - gStartPos - 1);
            }
            case None => "";
          }
        }
      }
    }
    return {(i : Int) => {
      if(i < offsetSeq.length) Some(offsetSeq.charAt(i))
      else None;
    }}
    
  }
  def saveToString() : String = {return geneID+"\t"+txID+"\t"+chrom+"\t"+strand+"\t"+gSpans.map{case (s,e) => s+","+e}.mkString(";") + "\t"+gCdsStart+"\t"+gCdsEnd+"\t"+seqGS.mkString("")+"\t"+leadingExonSeq.mkString(";")+"\t"+trailingExonSeq.mkString(";")+"\t"+gCdsNS._1+";"+gCdsNS._2;}
  
  def getRnaMut(ref : String, alt : String, gPos : Int) : String = {
    TXUtil.getShiftedMutString(ref,alt, pos = gPos, getPosString = {(g) => getRPos(g)}, swapStrand = (strand == '-'),
                 getFlankingSeq = getOffsetFlankSeqFunction(gPos)
                 );
  }
  def getCdsMut(ref : String, alt : String, gPos : Int) : String = {
    TXUtil.getShiftedMutString(ref,alt, pos = gPos, getPosString = {(g) => getCPos(g)}, swapStrand = (strand == '-'),
                 getFlankingSeq = getOffsetFlankSeqFunction(gPos)
                 );
  }
  
  def pVarInfo(pvar : String, start : Int, end : Int, alt : String, cType : String, sevType : String, pType : String, subType : String) : TXUtil.pVariantInfo = {
    TXUtil.pVariantInfo(this.txID,pvar,start,end,alt,cType,sevType,pType,subType);
  }

  def getProteinMut(ref : String, alt :String, gPos : Int, spliceRegionBuffer : Int = 2, recursionDepth : Int = 0) : TXUtil.pVariantInfo = {
    val (gVarStart,gVarEnd) = (gPos,gPos+ ref.length);
    
    val refS = if(strand == '+') ref else complementSequenceString(ref.reverse).mkString("");
    val altS = if(strand == '+') alt else complementSequenceString(alt.reverse).mkString("");
    val cPosOpt = convertGtoC(gPos);
    
    val vType = if(ref.length == 1 && alt.length == 1) "SNV"
                   else if(ref.length == 0) "INS";
                   else if(alt.length == 0) "DEL";
                   else                     "INDEL";
    
    
     
    if(ref == alt){
      //return ("p.=","NOVAR_SYNON_NOCHANGE");
      //return ("p.=",TXUtil.pVariantInfo(this.txID, -1, -1, "", "NOVAR","SYNON","NOCHANGE",""));
      pVarInfo("p.=",-1,-1,"","NOVAR","SYNON","NOCHANGE","NOCHANGE");
    } else if(ref.length > 0 && alt.length > 0 && ref.head == alt.head){
      val matchIdx = ref.zip(alt).indexWhere{case (r,a) => { r != a }}
      val matchIdxFinal = if(matchIdx == -1) math.min(ref.length,alt.length) else matchIdx
      if(recursionDepth > 10) reportln("   calling deep recursion ("+(recursionDepth+1)+"):\n"+
                                       "      getProteinMut(\""+ref+"\".drop("+matchIdxFinal+"), \""+alt+"\".drop("+matchIdxFinal+"),"+gPos+"+"+ matchIdxFinal+")\n"+
                                       "      getProteinMut(\""+ref.drop(matchIdxFinal)+"\",\""+alt.drop(matchIdxFinal)+"\","+(gPos + matchIdxFinal)+")",
                                       "warn");
      return getProteinMut(ref.drop(matchIdxFinal),alt.drop(matchIdxFinal),gPos + matchIdxFinal, recursionDepth = recursionDepth + 1);
    } else if(ref.length > 0 && alt.length > 0 && ref.last == alt.last){
      val matchIdx = ref.reverse.zip(alt.reverse).indexWhere{case (r,a) => { r != a }}
      val matchIdxFinal = if(matchIdx == -1) math.min(ref.length,alt.length) else matchIdx;
      if(recursionDepth > 10) reportln("   calling deep recursion:\n"+
                                       "      getProteinMut(\""+ref+"\".dropRight("+matchIdxFinal+"), \""+alt+"\".dropRight("+matchIdxFinal+"),"+(gPos)+")\n"+
                                       "      getProteinMut(\""+ref.dropRight(matchIdxFinal)+"\",\""+alt.dropRight(matchIdxFinal)+"\","+(gPos)+")","warn");
      return getProteinMut(ref.dropRight(matchIdxFinal),alt.dropRight(matchIdxFinal),gPos, recursionDepth = recursionDepth + 1);
    } else if(ref.length == 1 && alt.length == 1){
      if(gPos < gSpans.head._1){
        //if(strand == '+') return ("p.=",vType+"_PSYNON_nearGene-5");
        //else              return ("p.=",vType+"_PSYNON_nearGene-3");
        if(strand == '+') return pVarInfo("p.=",0,0,                        "",vType,"PSYNON","nearGene-5","nearGene-5");
        else              return pVarInfo("p.=",codons.length,codons.length,"",vType,"PSYNON","nearGene-3","nearGene-3");
      } else if(gPos >= gSpans.last._2){
        //if(strand == '+') return ("p.=",vType+"_PSYNON_nearGene-3");
        //else              return ("p.=",vType+"_PSYNON_nearGene-5");
        if(strand == '+') return pVarInfo("p.=",codons.length,codons.length,"",vType,"PSYNON","nearGene-3","nearGene-3");
        else              return pVarInfo("p.=",0,0,                        "",vType,"PSYNON","nearGene-5","nearGene-5");
      } else if(  (strand == '+' && gPos < gCdsStart) || 
                  (strand == '-' && gPos >= gCdsEnd)){
        //return ("p.?[UTR5p]",vType+"_UNK_UTR-5");
        return pVarInfo("p.?[UTR5p]",0,0,                        "",vType,"UNK","UTR-5","");
      } else if(  (strand == '+' && gPos >= gCdsEnd) || 
                  (strand == '-' && gPos <  gCdsStart)){
        //return ("p.?[UTR3p]",vType+"_UNK_UTR-3");
        return pVarInfo("p.?[UTR3p]",codons.length,codons.length,"",vType,"UNK","UTR-3","");
      } else if(cPosOpt.isEmpty && gSpans.exists{case (s,e) => { gPos < s && gPos + spliceRegionBuffer >= s }}){
        val lastUnchangedIdx = Range(0,cLen).find{ i => { val g = convertCtoG(i).get; (strand == '+' && g >= gPos) || (strand == '-' && g < gPos) }}.getOrElse(-3) / 3
        //if( strand == '+') return ("p.?[snvPossBreakSJDonor"+(lastUnchangedIdx+1)+"]",vType+"_PLOF_splice-3_PossibleSplice3Loss");
        //else               return ("p.?[snvPossBreakSJAccpt"+(lastUnchangedIdx+1)+"]",vType+"_PLOF_splice-5_PossibleSplice5Loss");
        if( strand == '+')  return pVarInfo("p.?[snvPossBreakSJDonor"+(lastUnchangedIdx+1)+"]",lastUnchangedIdx,lastUnchangedIdx,"",vType,"PLOF","splice-3","PossibleSplice3Loss");
        else                return pVarInfo("p.?[snvPossBreakSJDonor"+(lastUnchangedIdx+1)+"]",lastUnchangedIdx,lastUnchangedIdx,"",vType,"PLOF","splice-5","PossibleSplice5Loss");
      } else if(cPosOpt.isEmpty && gSpans.exists{case (s,e) => { gPos >= e && gPos - spliceRegionBuffer < e }}){
        val lastUnchangedIdx = Range(0,cLen).find{ i => { val g = convertCtoG(i).get; (strand == '+' && g >= gPos) || (strand == '-' && g < gPos) }}.getOrElse(-3) / 3
        if( strand == '-')  return pVarInfo("p.?[snvPossBreakSJDonor"+(lastUnchangedIdx+1)+"]",lastUnchangedIdx,lastUnchangedIdx,"",vType,"PLOF","splice-3","PossibleSplice3Loss");
        else                return pVarInfo("p.?[snvPossBreakSJDonor"+(lastUnchangedIdx+1)+"]",lastUnchangedIdx,lastUnchangedIdx,"",vType,"PLOF","splice-5","PossibleSplice5Loss");
        
      } else if(cPosOpt.isEmpty && gIntronSpans.exists{case (s,e) => { gPos < e && gPos >= s }}){
        val lastUnchangedIdx = Range(0,cLen).find{ i => { val g = convertCtoG(i).get; (strand == '+' && g >= gPos) || (strand == '-' && g < gPos) }}.getOrElse(-3) / 3;
        //return ("p.?[intronic"+(lastUnchangedIdx+1)+"]",vType+"_PSYNON_intron");
        return pVarInfo("p.?[intronic"+(lastUnchangedIdx+1)+"]",lastUnchangedIdx,lastUnchangedIdx,"",vType,"PSYNON","intron","intron");
      }
      
      convertGtoC(gPos) match {
        case Some(c) => {
          val codonIdx = c/3;
          val refCodon = codons(codonIdx);
          val refAA = getAminoAcidFromCodon(refCodon,aaNameStyle=1);
          val altCodon = TXUtil.swapCodonNS(refCodon,altS,c)
          val altAA = getAminoAcidFromCodon(altCodon,aaNameStyle=1);
          
          if(refAA == altAA){
            //return ("p."+refAA+(codonIdx+1)+"=",vType+"_SYNON_cds-synon");
            return pVarInfo("p."+refAA+(codonIdx+1)+"=",codonIdx,codonIdx+1,altAA,vType,"SYNON","cds-synon","cds-synon");
          } else if(refAA != "*" && altAA == "*"){
            //return ("p."+refAA+(codonIdx+1)+"*",vType+"_LLOF_STOP-GAIN");
            return pVarInfo("p."+refAA+(codonIdx+1)+"*",codonIdx,codonIdx+1,"*",vType,"LLOF","STOP-GAIN","STOP-GAIN");
          } else if(refAA == "*" && altAA != "*"){
            //3' Elongation!
            //return ("p."+refAA+(codonIdx+1)+altAA+"ext*?",vType+"_NONSYNON_STOP-LOSS");
            return pVarInfo("p."+refAA+(codonIdx+1)+altAA+"ext*?",codonIdx,codonIdx+1,altAA,vType,"NONSYNON","STOP-LOSS","STOP-LOSS");
          } else if(codonIdx == 0 && refAA == "Met" && altAA != "Met"){
            //NONSTANDARD:
            //return ("p."+refAA+(codonIdx+1)+altAA+"?",vType+"_PLOF_START-LOSS");
            return pVarInfo("p."+refAA+(codonIdx+1)+altAA+"?",codonIdx,codonIdx+1,altAA,vType,"PLOF","START-LOSS","START-LOSS");
          } else {
            //return ("p."+refAA+(codonIdx+1)+altAA,vType+"_NONSYNON_missense_swapAA");
            return pVarInfo("p."+refAA+(codonIdx+1)+altAA,codonIdx,codonIdx+1,altAA,vType,"NONSYNON","missense","swapAA");
          }
        }
        case None => {
          //return ("p.???",vType+"_???_???");
          return pVarInfo("p.???",-1,-1,"",vType,"UNK","???","???");
        }
      }
    } else {
      if(gVarStart <= gSpans.head._1 && gSpans.last._2 < gVarEnd){
        //ins/del covers entire gene:
        //return ("p.0",vType+"_LLOF_total-loss");
        return  pVarInfo("p.0",-1,codons.length,"",vType,"LLOF","total-loss","total-loss");
      } else if(gVarEnd <= gSpans.head._1){
        //if(strand == '+') return ("p.=",vType+"_PSYNON_nearGene-5");
        //else              return ("p.=",vType+"_PSYNON_nearGene-3");
        if(strand == '+') return pVarInfo("p.=",-1,-1,                      "",vType,"PSYNON","nearGene-5","nearGene");
        else              return pVarInfo("p.=",codons.length,codons.length,"",vType,"PSYNON","nearGene-3","nearGene");
      } else if(gVarStart >= gSpans.last._2){
        //if(strand == '+') return ("p.=",vType+"_PSYNON_nearGene-3");
        //else              return ("p.=",vType+"_PSYNON_nearGene-5");
        if(strand == '-') return pVarInfo("p.=",-1,-1,                      "",vType,"PSYNON","nearGene-5","nearGene");
        else              return pVarInfo("p.=",codons.length,codons.length,"",vType,"PSYNON","nearGene-3","nearGene");
      } else if( (strand == '+' && gVarEnd <= gCdsStart) || 
                 (strand == '-' && gVarStart >= gCdsEnd)){
          //return ("p.?[UTR5p]",vType+"_UNK_UTR-5");
          return pVarInfo("p.?[UTR5p]",-1,-1,                      "",vType,"UNK","UTR-5","UTR");
      } else if( (strand == '+' && (gVarStart >= gCdsEnd)) || 
                 (strand == '-' && (gVarEnd <= gCdsStart))){
          //return ("p.?[UTR3p]",vType+"_UNK_UTR-3");
        return pVarInfo("p.?[UTR3p]",codons.length,codons.length,"",vType,"UNK","UTR-3","UTR");
      } else {
        val lastUnchangedIdxOpt = Range(0,cLen).find{ i => {
          val g = convertCtoG(i).get
          (strand == '+' && g >= gVarEnd) || (strand == '-' && g < gVarStart)
        }}
        val lastUnchangedCodonNum = lastUnchangedIdxOpt match {
          case Some(idx) => idx/3 + 1;
          case None => cLen/3+1;
        }
        
        if(gSpansCDS.exists{case (s,e) => { gVarStart >= s && gVarEnd <= e && gVarStart < e && gVarEnd > s }}){
          //in/del spans a single coding-region (exon):
          val cVarSpan = Vector(convertGtoC(gVarStart).get,convertGtoC(math.max(gVarStart,gVarEnd-1)).get)
          val (cVarStart,cVarLast) = (cVarSpan.min, cVarSpan.max);
          val (aVarStart,aVarLast) = (cVarStart/3,cVarLast/3);
          
          val codonStartOffset = cVarStart % 3
          
          val postVarSeq = if(ref.length == 0) cSeq.drop(cVarLast) else cSeq.drop(cVarLast+1)
          
          val refCodons = codons.drop(aVarStart) //Range(aVarStart,codons.length).map(codons(_));
          val altCodons = getCodonsFromSeq( (codons(aVarStart).substring(0,codonStartOffset) + altS).toVector ++ postVarSeq )
          val refAA = refCodons.map(getAminoAcidFromCodon(_,aaNameStyle=1)).zip(Stream.from(aVarStart));
          val altAA = altCodons.map(getAminoAcidFromCodon(_,aaNameStyle=1));
          
          val leadZip = refAA.zip(altAA)
          val leadMatch = leadZip.takeWhile{case ((r,i),a) => {r == a}}
          
          val lastMatchingIdx = aVarStart + leadMatch.length - 1;
          val firstMismatchingIdx = aVarStart + leadMatch.length;
          
          if(leadMatch.length == refAA.length){
            //return ("p.=",vType+"_SYNON_cds-synon_SYNONINDEL");
            return pVarInfo("p."+refAA.head._1+aVarStart+"=",aVarStart,aVarStart+1,altAA.head,vType,"SYNON","cds-synon","cds-synon");
          }
          if(leadMatch.length == refAA.length - 1){
            val altCodonsExt = getCodonsFromSeq((codons(aVarStart).substring(0,codonStartOffset) + altS).toVector ++ postVarSeq ++ cTrailingSeq.toVector);
            val altAAext = altCodonsExt.map(getAminoAcidFromCodon(_,aaNameStyle=1));
            //if(altAAext.contains("*")){
            //  return ("p.*"+(firstMismatchingIdx+1)+altAAext(firstMismatchingIdx)+"ext"+,"INDEL_STOPLOSS");
            //} else {
              //return ("p.*"+(firstMismatchingIdx+1)+altAAext(leadMatch.length)+"ext*?",vType+"_NONSYNON_STOP-LOSS");
            return pVarInfo("p.*"+(firstMismatchingIdx+1)+altAAext(leadMatch.length)+"ext*?",firstMismatchingIdx,firstMismatchingIdx+1,altAAext(leadMatch.length),vType,"NONSYNON","STOP-LOSS","STOP-LOSS");
            //}
          }
          if(refAA.length > 0 && refAA.head._2 == 0 && leadMatch.length == 0){
            //Start loss?
            //return ("p.Met1"+altAA.head+"?",vType+"_PLOF_START-LOSS");
            return pVarInfo("p.Met1"+altAA.head+"?",0,1,altAA.head,vType,"PLOF","START-LOSS","START-LOSS");
          }
          
          if(math.abs(ref.length - alt.length) % 3 != 0){
            //Frameshift:
            val altAAext = getCodonsFromSeq((codons(aVarStart).substring(0,codonStartOffset) + altS).toVector ++ postVarSeq ++ cTrailingSeq).map(getAminoAcidFromCodon(_,aaNameStyle=1));
            if(aminos(firstMismatchingIdx) == "*"){
              //return ("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+altAAext(firstMismatchingIdx-aVarStart)+"ext*?",vType+"_NONSYNON_STOP-LOSS");
              return pVarInfo("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+altAAext(firstMismatchingIdx-aVarStart)+"ext*?",firstMismatchingIdx,firstMismatchingIdx+1,
                              altAAext(firstMismatchingIdx-aVarStart),vType,"NONSYNON","STOP-LOSS","STOP-LOSS");
            } else if(altAAext(firstMismatchingIdx-aVarStart) == "*"){
              //return ("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+"*",vType+"_LLOF_STOP-GAIN_fsSTOP");
              return pVarInfo("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+"*",firstMismatchingIdx,firstMismatchingIdx+1,
                              "*",vType,"LLOF","STOP-GAIN","fsSTOP");
            } else {
              //return ("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+altAAext(firstMismatchingIdx-aVarStart)+"fs",vType+"_LLOF_frameshift");
              return pVarInfo("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+altAAext(firstMismatchingIdx-aVarStart)+"fs",firstMismatchingIdx,firstMismatchingIdx+1,
                              altAAext(firstMismatchingIdx-aVarStart),vType,"LLOF","frameshift","fs");
            }
          }
          
          val tailMatch = refAA.drop(leadMatch.length).reverse.zip(altAA.drop(leadMatch.length).reverse).takeWhile{case ((r,i),a) => r == a}
          //val tailMatchingIdx = refAA.last._2 - tailMatch.length + 1;
          //val tailMismatchingIdx = refAA.last._2 - tailMatch.length;
          
          val mmRef = refAA.slice(leadMatch.length,refAA.length - tailMatch.length);
          val mmAlt = altAA.slice(leadMatch.length,altAA.length - tailMatch.length);
          
          if(mmAlt.contains("*")){
            val stopCodonIdx = mmAlt.indexOf("*");
            val mmAltSub = mmAlt.take(stopCodonIdx+1);
            if(mmAltSub.length == 1){
               //return ("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+"*",vType+"_LLOF_STOP-GAIN");
               return pVarInfo("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+"*",firstMismatchingIdx,firstMismatchingIdx+1,
                               "*",vType,"LLOF","STOP-GAIN","STOP-GAIN");
            } else {
              //return ("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+"delins"+mmAltSub.mkString(""),vType+"_LLOF_STOP-GAIN_delinsToStopGain");
               return pVarInfo("p."+aminos(firstMismatchingIdx)+(firstMismatchingIdx+1)+"delins"+mmAltSub.mkString(""),firstMismatchingIdx,firstMismatchingIdx+1,
                               mmAltSub.mkString(""),vType,"LLOF","STOP-GAIN","delinsToStopGain");
            }
          } else if(mmAlt.length == 0){
            //deletion:
            if(mmRef.length == 1){
              //return ("p."+mmRef.head._1+(mmRef.head._2+1)+"del",vType+"_NONSYNON_missense_delAA")
               return pVarInfo("p."+mmRef.head._1+(mmRef.head._2+1)+"del",mmRef.head._2,mmRef.head._2+1,
                               "",vType,"NONSYNON","missense","delAA");
            } else {
              //return ("p."+mmRef.head._1+(mmRef.head._2+1)+"_"+mmRef.last._1+(mmRef.last._2+1)+"del",vType+"_NONSYNON_missense_delAA")
               return pVarInfo("p."+mmRef.head._1+(mmRef.head._2+1)+"del",mmRef.head._2,mmRef.last._2+1,
                               "",vType,"NONSYNON","missense","delAA");
            }
          } else if(mmRef.length == 0){
            //insertion:
            //return ("p."+aminos(lastMatchingIdx)+(lastMatchingIdx+1)+"_"+aminos(lastMatchingIdx+1)+(lastMatchingIdx+2)+"ins"+mmAlt.mkString(""),vType+"_NONSYNON_missense_insAA")
             return pVarInfo("p."+aminos(lastMatchingIdx)+(lastMatchingIdx+1)+"_"+aminos(lastMatchingIdx+1)+(lastMatchingIdx+2)+"ins"+mmAlt.mkString(""),lastMatchingIdx,lastMatchingIdx+1,
                             mmAlt.mkString(""),vType,"NONSYNON","missense","insAA");
          } else if(mmRef.length == 1 && mmAlt.length == 1){
            //return ("p."+mmRef.head._1+(mmRef.head._2+1)+mmAlt.head,vType+"_NONSYNON_missense_swapAA")
             return pVarInfo("p."+mmRef.head._1+(mmRef.head._2+1)+mmAlt.head,mmRef.head._2,mmRef.head._2+1,
                             mmAlt.head,vType,"NONSYNON","missense","swapAA");
          } else if(mmRef.length == 1){
            //single-codon delins:
            //return ("p."+mmRef.head._1+(mmRef.head._2+1)+"delins"+mmAlt.mkString(""),vType+"_NONSYNON_missense_indelAA")
             return pVarInfo("p."+mmRef.head._1+(mmRef.head._2+1)+"delins"+mmAlt.mkString(""),mmRef.head._2,mmRef.head._2+1,
                             mmAlt.mkString(""),vType,"NONSYNON","missense","indelAA");
          } else {
            //multi-codon delins:
            //if(mmRef.length == mmAlt.length) return ("p."+mmRef.head._1+(mmRef.head._2+1)+"_"+mmRef.last._1+(mmRef.last._2+1)+"delins"+mmAlt.mkString(""),vType+"_NONSYNON_missense_multSwapAA");
            //else                             return ("p."+mmRef.head._1+(mmRef.head._2+1)+"_"+mmRef.last._1+(mmRef.last._2+1)+"delins"+mmAlt.mkString(""),vType+"_NONSYNON_missense_indelAA");
            if(mmRef.length == mmAlt.length){
              return pVarInfo("p."+mmRef.head._1+(mmRef.head._2+1)+"_"+mmRef.last._1+(mmRef.last._2+1)+"delins"+mmAlt.mkString(""),mmRef.head._2,mmRef.last._2+1,
                               mmAlt.mkString(""),vType,"NONSYNON","missense","multSwapAA");
            } else {
              return pVarInfo("p."+mmRef.head._1+(mmRef.head._2+1)+"_"+mmRef.last._1+(mmRef.last._2+1)+"delins"+mmAlt.mkString(""),mmRef.head._2,mmRef.last._2+1,
                               mmAlt.mkString(""),vType,"NONSYNON","missense","indelAA");
            }

            
          }
        } else if(gSpans.exists{case (s,e) => { gVarStart > s && gVarStart < e}} && gSpans.exists{case (s,e) => { gVarEnd > s && gVarEnd < e}}){
          //in/del starts and ends on different exons:
          //return ("p.?[FullIntronInsOrDel"+(lastUnchangedCodonNum+1)+"]",vType+"_UNK_FullIntronIndel");
          return pVarInfo("p.?[FullIntronInsOrDel"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"UNK","FullIntronIndel","FullIntronIndel");
        } else if(gSpans.exists{case (s,e) => { gVarStart <= s && gVarEnd >= e}}){
          //in/del deletes or replaces an entire exon
          //return ("p.?[FullExonInsOrDel"+(lastUnchangedCodonNum+1)+"]",vType+"_PLOF_FullExonIndel");
          return pVarInfo("p.?[FullExonInsOrDel"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","FullExonIndel","FullExonIndel");
        } else if(gVarStart < gCdsStart && gVarEnd > gCdsStart){
          //if(strand == '+') return("p.?[IndelStartSite]",vType+"_PLOF_START-LOSS_startIndel");
          //else              return("p.?[IndelStopSite]", vType+"_NONSYNON_STOP-LOSS_stopIndel");
          if(strand == '+') return pVarInfo("p.?[IndelStartSite]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","START-LOSS","startIndel");
          else              return pVarInfo("p.?[IndelStopSite]", lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"NONSYNON","STOP-LOSS","stopIndel");
        } else if(gVarStart < gCdsEnd && gVarEnd > gCdsEnd){
          //if(strand == '+') return("p.?[IndelStopSite]", vType+"_NONSYNON_STOP-LOSS_stopIndel");
          //else              return("p.?[IndelStartSite]",vType+"_PLOF_START-LOSS_startIndel");
          if(strand == '-') return pVarInfo("p.?[IndelStartSite]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","START-LOSS","startIndel");
          else              return pVarInfo("p.?[IndelStopSite]", lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"NONSYNON","STOP-LOSS","stopIndel");
        } else if(gVarEnd <= gCdsStart){
          //if(strand == '+') return("p.?[Indel5primeEnd]",vType+"_UNK_UTR-5");
          //else              return("p.?[Indel3primeEnd]",vType+"_UNK_UTR-3");
          if(strand == '+') return pVarInfo("p.?[Indel5primeEnd]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"UNK","UTR-5","UTR");
          else              return pVarInfo("p.?[Indel3primeEnd]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"UNK","UTR-3","UTR");
        } else if(gVarStart >= gCdsEnd){
          //if(strand == '+') return("p.?[Indel3primeEnd]",vType+"_UNK_UTR-3");
          //else              return("p.?[Indel5primeEnd]",vType+"_UNK_UTR-5");
          if(strand == '-') return pVarInfo("p.?[Indel5primeEnd]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"UNK","UTR-5","UTR");
          else              return pVarInfo("p.?[Indel3primeEnd]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"UNK","UTR-3","UTR");
        } else if(gSpans.exists{case (s,e) => {gVarStart < e && gVarEnd > e}}){
          //if(strand == '+') return ("p.?[indelBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]", vType+"_PLOF_splice-5_Splice5Loss");
          //else              return ("p.?[indelBreakSJDonor"+(lastUnchangedCodonNum+1)+"]", vType+"_PLOF_splice-3_Splice3Loss");
          if(strand == '+') return pVarInfo("p.?[indelBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-5","splice");
          else              return pVarInfo("p.?[indelBreakSJDonor"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-3","splice");
        } else if(gSpans.exists{case (s,e) => {gVarStart < s && gVarEnd > s}}){
          //if(strand == '+') return ("p.?[indelBreakSJDonor"+(lastUnchangedCodonNum+1)+"]", vType+"_PLOF_splice-3_Splice3Loss");
          //else              return ("p.?[indelBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]", vType+"_PLOF_splice-5_Splice5Loss");
          if(strand == '-') return pVarInfo("p.?[indelBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-5","splice");
          else              return pVarInfo("p.?[indelBreakSJDonor"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-3","splice");
        } else if(gSpans.exists{case (s,e) => {gVarStart < e + spliceRegionBuffer && gVarEnd > e}}){
          //if(strand == '+') return ("p.?[indelPossBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]",  vType+"_PLOF_splice-5_PossibleSplice5Loss");
          //else              return ("p.?[indelPossBreakSJDonor"+(lastUnchangedCodonNum+1)+"]",  vType+"_PLOF_splice-3_PossibleSplice3Loss");
          if(strand == '+') return pVarInfo("p.?[indelPossBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-5","possSplice");
          else              return pVarInfo("p.?[indelPossBreakSJDonor"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-3","possSplice");
        } else if(gSpans.exists{case (s,e) => {gVarStart < s && gVarEnd > s - spliceRegionBuffer}}){
          //if(strand == '+') return ("p.?[indelPossBreakSJDonor"+(lastUnchangedCodonNum+1)+"]",  vType+"_PLOF_splice-3_PossibleSplice3Loss");
          //else              return ("p.?[indelPossBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]",  vType+"_PLOF_splice-5_PossibleSplice5Loss");
          if(strand == '-') return pVarInfo("p.?[indelPossBreakSJAccpt"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-5","possSplice");
          else              return pVarInfo("p.?[indelPossBreakSJDonor"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PLOF","splice-3","possSplice");
        } else if(gIntronSpans.exists{ case (s,e) => {s < gVarEnd && gVarStart < e}}){
          //return ("p.?[intronic"+(lastUnchangedCodonNum+1)+"]", vType+"_PSYNON_intron");
          return pVarInfo("p.?[intronic"+(lastUnchangedCodonNum+1)+"]",lastUnchangedCodonNum,lastUnchangedCodonNum,"",vType,"PSYNON","intron","intron");
        } else {
          //return ("p.???OtherVariantType",vType+"_UNK_UNKNOWN");
          return pVarInfo("p.???",-1,-1,"",vType,"UNK","???","???");
        }
      }
    }
    return pVarInfo("p.???",-1,-1,"",vType,"UNK","???","???");
  }
  
}

/**
 *  
 *  The output has the following format:
 *  (HGVS protein variant description, variantTypeCode)
 *  Where the HGVS protein variant description is defined as per the HGVS sequence variant nomenclature guidelines
 *  found here: (http://varnomen.hgvs.org/). Some variants may have additional information contained in square-brackets.
 *  
 *  The variantTypeCode is an underscore-delimited code consisting of three to four fields in the form:
 *  nvartype_severityCategory_consequence[_subtype]
 *  
 ** Nucleotide Variant types:
 *  SNV, INS, DEL, INDEL
 *  
 ** Amino Acid Effect Categories:
 *  LLOF: Likely-Loss-Of-Function variant. These are nonsense stop-gain variants,  
 *        frameshift variants, and large full-exon deletions.
 *  PLOF: Probable-Loss-Of-Function variant. These are splice site indels or loss 
 *        of the start codon.
 *  NONSYNON: Non-Synonymous. A variant that results in a change in the amino acid sequence.
 *            These include synonymous point mutations, simple in-frame indels, and 
 *            stop-loss variants.
 *  UNK: A variant that may or may not result in a change.
 *  PSYNON: Probably-Synonymous. A variant that probably has no effect, but could affect expression.
 *          Includes intronic variants and variants that are near genes.
 *  SYNON: Synonymous. A variant that results in no amino acid change.
 *  
 *  Molecular Consequence Types:
 *   splice-3
 *   splice-5
 *   STOP-LOSS
 *   missense
 *  STOP-GAIN
 *  frameshift
 *   UTR-5
 *   UTR-3
 *   intron
 *   nearGene-3
 *   nearGene-5
 *   cds-synon
 *  START-LOSS
 *  
 *  Molecular Consequence Subtypes (optional):
 *  PossibleSplice5Loss: This is present when the canonical 5prime splice site region is covered by a variant.
 *  PossibleSplice3Loss: This is present when the canonical 3prime splice site region is covered by a variant.
 *  Splice5Loss: This is present when an entire splice site is deleted or replaced.
 *  Splice3Loss: This is present when an entire splice site is deleted or replaced.
 *  insAA: This is a type of missense variant where one or more amino acids is inserted into the protein sequence.
 *  delAA: This is a type of missense variant where one or more amino acids is removed from the protein sequence.
 *  swapAA: This is a type of missense variant where exactly one amino acid is swapped for another.
 *  indelAA: This is a type of missense variant where one or more amino acids are exchanged for one or more different amino acids.
 *  
 */

  

  

/*
 * LOF types:
 * INDEL_SpliceAccptLoss,INDEL_SpliceDonorLoss,INDEL_PossibleSpliceAccptLoss,INDEL_PossibleSpliceDonorLoss,INDEL_StartSite,INDEL_StopSite,INDEL_DelinsToNonsense,INDEL_Nonsense
 */



object TXUtil {
  
  object pVariantInfo {
    sealed trait SeverityType;
    case object LLOF extends SeverityType;
    case object PLOF extends SeverityType;
    case object NONSYNON extends SeverityType;
    case object UNK extends SeverityType;
    case object PSYN extends SeverityType;
    case object SYNON extends SeverityType;
    
    sealed trait CType;
    case object SNV extends CType;
    case object INS extends CType;
    case object DEL extends CType;
    case object INDEL extends CType;
    
    sealed trait PType;
    case object splice3 extends PType;
    case object splice5 extends PType;
    case object STOPLOSS extends PType;
    case object missense extends PType;
    case object STOPGAIN extends PType;
    case object frameshift extends PType;
    case object UTR5 extends PType;
    case object UTR3 extends PType;
    case object intron extends PType;
    case object nearGene3 extends PType;
    case object neargene5 extends PType;
    case object cdssynon extends PType;
    case object STARTLOSS extends PType;
    
    //TO DO: switch this over to case class enums?
  }

  //pVariantInfo(txid, start, end, alt, varTypes : Seq[String])
  case class pVariantInfo(txid : String, pvar : String, start : Int, end : Int, altAA : String, cType : String, severityType : String, pType : String, subType : String, ID : String = ""){
    lazy val isSwap : Boolean = subType == "swapAA";
    lazy val varType : String = cType + "_" + severityType + "_" + pType +"_"+subType;
    def saveToString(delim : String = "/") : String = {
      txid + delim + pvar + delim + start + delim + end + delim + altAA + delim + cType + delim + severityType + delim + pType + delim + subType;
    }
  }
  def getPvarInfoFromString(str : String, delim : String = "/", ID : String = "") : pVariantInfo = {
    val cells = str.split(delim,-1);
    if(cells.length != 9){
      error("Error: pVariantInfo unpack: supplied string is malformatted: length = "+cells.length+"\n Supplied String=\""+str+"\"\ncells = [\""+cells.mkString("\",\"")+"\"]");
    }
    return pVariantInfo(txid = cells(0), 
                 pvar = cells(1), 
                 start = string2int(cells(2)), 
                 end  = string2int(cells(3)), 
                 altAA = cells(4), 
                 cType = cells(5), 
                 severityType = cells(6), 
                 pType = cells(7), 
                 subType = cells(8),
                 ID = ID);
    
  }
  
  
  val VARIANT_TYPE_LIST_IN_ORDER_TOP = Vector[String]("SYNON","PSYNON","UNK","NONSYNON","PLOF","LLOF");
  
    
  def getWorstProteinMut(varP : Seq[(String,String)], txList : Seq[String], topVarOrd : Vector[String] = VARIANT_TYPE_LIST_IN_ORDER_TOP) : (String,String,String) = {
    if(varP.length == 0) return ("","","");
    
    val (v,t,i) = varP.zip(txList).foldLeft(("",List[String](),-1)){case ((soFarVar,soFarTypeSet,soFarIdx),((currVar,currType),txID)) => {
      val currIdx = topVarOrd.indexOf( currType.split("_")(1) );
      if(soFarIdx < currIdx){
        (txID+":"+currVar,List[String](currType),currIdx);
      } else if(soFarIdx == currIdx){
        if(! soFarTypeSet.contains(currType)){
          (soFarVar + "|" + txID + ":" + currVar,soFarTypeSet :+ currType,soFarIdx);
        } else {
          (soFarVar,soFarTypeSet,soFarIdx);
        }
      } else {
        (soFarVar,soFarTypeSet,soFarIdx);
      }
    }}
    (v,t.sorted.mkString("|"),topVarOrd(i));
  }
  
  def swapCodon(codon : String, alt : String, cPos : Int, strand : Char) : String = {
    def ALT = if(strand == '+') alt else complementSequenceString(alt.toSeq).head.toString();
    if(cPos % 3 == 0){
      return ALT + codon.substring(1,3);
    } else if(cPos % 3 == 1){
      return codon.head + ALT + codon.last;
    } else {
      return codon.init + ALT;
    }
  }
  def swapCodonNS(codon : String, altS : String, cPos : Int) : String = {
    if(cPos % 3 == 0){
      return altS + codon.substring(1,3);
    } else if(cPos % 3 == 1){
      return codon.head + altS + codon.last;
    } else {
      return codon.init + altS;
    }
  }
  //def delCodon(cSeq : Vector[Char], len : Int, cPos : Int) : String = {
    
  //}

  case class TranscriptHolder(geneID : String, txID : String, chrom : String, strand : Char, gSpans : Vector[(Int,Int)], 
                              gCdsStart : Int = Int.MaxValue, gCdsEnd : Int = Int.MinValue, 
                              exonSeqGS : Vector[Vector[Char]],
                              leadingExonSeq : Vector[String],
                              trailingExonSeq : Vector[String],
                              bufferSize : Int = 32){
    
    def complete(addStopCodon : Boolean = true) : Option[TXUtil] = {
      if(! isValid){
        if(gSpans.isEmpty) warning("Invalid TX (No Exons): " + txID,"warning.InvalidTX.NoExon",100);
        else warning("Invalid TX (No CDS): " + txID,"warning.InvalidTX.NoCDS",100);
        return None;
        //error("Invalid TX! " + txID);
      }
      var gSpansOut = gSpans;
      var gCdsStartOut = gCdsStart;
      var gCdsEndOut = gCdsEnd;
      var exonSeq = exonSeqGS.flatten;
      var leadingExonSeqOut = leadingExonSeq.map(s => {repString("N",bufferSize - s.length) + s});
      var trailingExonSeqOut = trailingExonSeq.map(s => {s + repString("N",bufferSize - s.length)});

      if(addStopCodon){
        if(strand == '+'){
          if(gSpans.last._2 < gCdsEndOut + 3){
            gCdsEndOut = gCdsEndOut + 3;
            val diff = gCdsEndOut - gSpansOut.last._2;
            val lastTrail = trailingExonSeqOut.last
            exonSeq = exonSeq ++ lastTrail.take(diff).toVector;
            trailingExonSeqOut = trailingExonSeqOut.updated(trailingExonSeqOut.length-1, lastTrail.drop(diff) + repString("N",diff));
            gSpansOut = gSpansOut.init :+ ((gSpansOut.last._1,gCdsEndOut));
          } else if(gSpans.exists{case (s,e) => {gCdsEndOut + 3 > e && gCdsEndOut <= e}}){
            val exonidx = gSpans.indexWhere{case (s,e) => {gCdsEndOut + 3 > e && gCdsEndOut <= e}}
            val exonend = gSpans(exonidx)._2
            if(exonend - gCdsEndOut == 0){
              gCdsEndOut = gSpans(exonidx+1)._1 + 3
            } else if(exonend - gCdsEndOut == 1){
              gCdsEndOut = gSpans(exonidx+1)._1 + 2
            } else {
              gCdsEndOut = gSpans(exonidx+1)._1 + 1
            }
          } else {
            gCdsEndOut = gCdsEndOut + 3;
          }
        } else {
          if(gSpans.head._1 > gCdsStartOut - 3){
            gCdsStartOut = gCdsStartOut - 3;
            val diff =  gSpans.head._1 - gCdsStartOut;
            val lastTrail = leadingExonSeqOut.head;
            exonSeq = lastTrail.takeRight(diff).toVector ++ exonSeq;
            leadingExonSeqOut = leadingExonSeqOut.updated(0, repString("N",diff) + lastTrail.dropRight(diff));
            gSpansOut = ((gCdsStartOut,gSpansOut.head._2)) +: gSpansOut.tail;
          } else if(gSpans.exists{case (s,e) => {gCdsStartOut - 3 < s && gCdsStartOut >= s}}){
            val exonidx = gSpans.indexWhere{case (s,e) => {gCdsStartOut - 3 < s && gCdsStartOut >= s}}
            val exonend = gSpans(exonidx)._1
            if(gCdsStartOut - exonend == 0){
              gCdsStartOut = gSpans(exonidx-1)._2 - 3
            } else if(gCdsStartOut - exonend == 1){
              gCdsStartOut = gSpans(exonidx-1)._2 - 2
            } else {
              gCdsStartOut = gSpans(exonidx-1)._2 - 1
            }
          } else {
            gCdsStartOut = gCdsStartOut - 3;
          }
        }
      }
      
      return Some(new TXUtil(geneID,txID,chrom,strand,
                        gSpansOut,
                        gCdsStartOut,
                        gCdsEndOut,
                        exonSeq, 
                        leadingExonSeqOut, 
                        trailingExonSeq,
                        gCdsNS = (gCdsStart,gCdsEnd)
                        ));
    }
    
    def isValid : Boolean = {
      gCdsStart != Int.MaxValue && gCdsEnd != Int.MinValue && (! gSpans.isEmpty)
    }
    def addSpan(start : Int, end : Int, seq : internalUtils.genomicAnnoUtils.EfficientGenomeSeqContainer) : TranscriptHolder = {
      val (l,t) = getBufferSeq(chrom=chrom,start,end,seq,bufferSize);
      val currSeq = seq.getSeqForInterval(chrom,start,end).toVector

      if(gSpans.length == 0){
        return TranscriptHolder(geneID, txID, chrom, strand,gSpans = Vector[(Int,Int)]((start,end)), 
                              gCdsStart=gCdsStart, gCdsEnd=gCdsEnd, exonSeqGS = Vector[Vector[Char]](currSeq),
                              leadingExonSeq = Vector(l),
                              trailingExonSeq = Vector(t),
                              bufferSize=bufferSize
                              );
      }
      
      if(  gSpans.exists{case (s,e) => { s < end && start < e }}  ){
        error("ERROR: Transcript Holder: BAD SPAN:\n"+
                           "   Attempted span addition of ("+start+","+end+")\n"+
                           "   For TX: "+gSpans.map{case (i,j) => "("+i+","+j+")"}.mkString(" "));
      }
      
      val (newLes,newTes, newSeq) = if(end <= gSpans.head._1){
        (l +: leadingExonSeq,t +: trailingExonSeq,currSeq +: exonSeqGS);
      } else if(start >= gSpans.last._2){
        (leadingExonSeq :+ l, trailingExonSeq :+ t, exonSeqGS :+ currSeq);
      } else {
        val exonIdx = gSpans.indexWhere{case (s,e) => { s > end }}
        if(exonIdx == -1) error("ERROR: Transcript Holder: BAD SPAN: Impossible state!");
        (leadingExonSeq.take(exonIdx) ++ Vector(l) ++ leadingExonSeq.drop(exonIdx),
         trailingExonSeq.take(exonIdx) ++ Vector(t) ++ trailingExonSeq.drop(exonIdx),
         exonSeqGS.take(exonIdx) ++ Vector[Vector[Char]](currSeq) ++ exonSeqGS.drop(exonIdx)
         )
      }
      
      return TranscriptHolder(geneID, txID, chrom, strand, gSpans :+ (start,end), 
                              gCdsStart=gCdsStart, gCdsEnd=gCdsEnd, exonSeqGS = newSeq,
                              leadingExonSeq = newLes,
                              trailingExonSeq = newTes,
                              bufferSize=bufferSize
                              );
    }
    def addCds(start : Int, end : Int) : TranscriptHolder = {
      TranscriptHolder(geneID, txID, chrom, strand, gSpans, gCdsStart=math.min(start,gCdsStart), gCdsEnd=math.max(end,gCdsEnd), exonSeqGS = exonSeqGS,
                      leadingExonSeq = leadingExonSeq, trailingExonSeq=trailingExonSeq,bufferSize = bufferSize);
    }
    def addCdsStart(x : Int) : TranscriptHolder = {
      TranscriptHolder(geneID, txID, chrom, strand, gSpans, gCdsStart=math.min(x,gCdsStart), gCdsEnd=gCdsEnd, exonSeqGS = exonSeqGS,
                       leadingExonSeq = leadingExonSeq, trailingExonSeq=trailingExonSeq,bufferSize = bufferSize);
    }
    def addCdsEnd(x : Int) : TranscriptHolder = {
      TranscriptHolder(geneID, txID, chrom, strand, gSpans, gCdsStart, gCdsEnd=math.max(x,gCdsEnd), exonSeqGS = exonSeqGS,
                       leadingExonSeq = leadingExonSeq, trailingExonSeq=trailingExonSeq,bufferSize = bufferSize);
    }
  }
  
  def getBufferSeq(chrom : String, start : Int, end : Int, genomeSeq : internalUtils.genomicAnnoUtils.EfficientGenomeSeqContainer, bufferSize : Int = 32) : (String,String) = {
    val bs = math.max(start-bufferSize,0)
    val be = math.min(end+bufferSize,genomeSeq.extendBufferToOrTruncate(end+bufferSize));
    
    val startN = bs - (start - bufferSize);
    val endN = (end + bufferSize) - be;
    
    ( repString("N",startN) + genomeSeq.getSeqForInterval(chrom,bs,start),
      genomeSeq.getSeqForInterval(chrom,end,be) + repString("N",endN) );
  }
  
  def createTranscriptHolder(gtfLine : StdGtfLine, genomeSeq : internalUtils.genomicAnnoUtils.EfficientGenomeSeqContainer, bufferSize : Int = 32) : TranscriptHolder = {
    if(gtfLine.isExon){
      val (start,end) = (gtfLine.start - 1,gtfLine.end)
      val (les,tes) = getBufferSeq(chrom=gtfLine.chromName,start,end,genomeSeq,bufferSize);
      
      TranscriptHolder(geneID = gtfLine.getGeneID, 
                       txID = gtfLine.getTxID, 
                       chrom = gtfLine.chromName, 
                       strand = gtfLine.strand, 
                       gSpans = Vector((start,end)), 
                       exonSeqGS = Vector[Vector[Char]](genomeSeq.getSeqForInterval(gtfLine.chromName,start,end).toVector),
                       leadingExonSeq = Vector(les),
                       trailingExonSeq = Vector(tes),
                       bufferSize = bufferSize
                       )
    } else {
      val (start,end) = (gtfLine.start - 1,gtfLine.end)
      TranscriptHolder(geneID = gtfLine.getGeneID, 
                       txID = gtfLine.getTxID, 
                       chrom = gtfLine.chromName, 
                       strand = gtfLine.strand, 
                       gSpans = Vector[(Int,Int)](), 
                       gCdsStart = start, gCdsEnd = end,
                       exonSeqGS = Vector[Vector[Char]](),
                       leadingExonSeq = Vector[String](),
                       trailingExonSeq = Vector[String](),
                       bufferSize = bufferSize
                       )
    }
  }
  
  var debugCount = 0;
  
  def buildTXUtilsFromAnnotation(gtffile : String, genomeFa : String, addStopCodon : Boolean = true, debugMode : Boolean = false, chromSet : Option[Set[String]] = None) : scala.collection.mutable.Map[String,TXUtil] = {
    val stdGtfCodes = new internalUtils.GtfTool.GtfCodes();
    val gtfReader = GtfReader.getStdGtfReader(gtffile, true, true, "\\s+", stdGtfCodes);
    
    val ipr = internalUtils.stdUtils.IteratorProgressReporter_ThreeLevel("lines", 2000, 10000 , 20000 )
    val wrappedGtfIterator = internalUtils.stdUtils.wrapIteratorWithProgressReporter(gtfReader , ipr )
    
    val genomeSeq = internalUtils.genomicAnnoUtils.buildEfficientGenomeSeqContainer(Seq(genomeFa));
     
    val holdermap : scala.collection.mutable.AnyRefMap[String,TranscriptHolder] = scala.collection.mutable.AnyRefMap[String,TranscriptHolder]();
    
    reportln("Starting TX parse:","progress")
    for(gtfLine <- wrappedGtfIterator){
      if((chromSet.isEmpty || chromSet.get.contains(gtfLine.chromName)) && (gtfLine.isCDS || gtfLine.isExon)){
      //val biotype = gtfLine.attributeMap.getOrElse("transcript_biotype","NA");
      //if(biotype == "protein_coding"){
      //if(gtfLine.chromName == "chr20" && debugCount < 50){
      //  genomeSeq.makeVerbose;
      //  debugCount = debugCount + 1;
      //} else {
      //  genomeSeq.makeQuiet;
      //}
        
      val tx = gtfLine.getTxID;
      genomeSeq.shiftBufferTo(gtfLine.chromName,gtfLine.start - 20000);
      
      holdermap.get( tx ) match {
        case Some(txHolder) => {
          if(gtfLine.isCDS){
            if(debugMode) reportln("CDS_ADDING"+"\t\""+tx+"\"\t"+gtfLine.str,"debug");
            holdermap += ((tx,txHolder.addCds(gtfLine.start-1,gtfLine.end)));
          } else if(gtfLine.isExon){
            if(debugMode) reportln("EXON_ADDING"+"\t\""+tx+"\"\t"+gtfLine.str,"debug");
            holdermap += ((tx,txHolder.addSpan(gtfLine.start-1,gtfLine.end,genomeSeq)));
          }
        }
        case None => {
          if(gtfLine.isCDS){
            //error("CDS cannot occur before first exon!");
            if(debugMode) reportln("CDS_NEW"+"\t\""+tx+"\"\t"+gtfLine.str,"debug");
            holdermap += ((tx, createTranscriptHolder(gtfLine, genomeSeq)));
          } else if(gtfLine.isExon){
            if(debugMode) reportln("EXON_NEW"+"\t\""+tx+"\"\t"+gtfLine.str,"debug");
            holdermap += ((tx, createTranscriptHolder(gtfLine, genomeSeq)));
          }
        }
      }} else {
        if(debugMode) reportln("SKIP"+"\t"+gtfLine.str,"debug");
      }
    }
    
    reportln("TX parse complete.","progress")
    
    reportln("Building TX Map.","progress")
    val out : scala.collection.mutable.AnyRefMap[String,TXUtil] = scala.collection.mutable.AnyRefMap[String,TXUtil]();
    for((tx,holder) <- holdermap.iterator){
      //if(holder.isValid){
      val txc = holder.complete(addStopCodon);
      if(txc.isDefined){
        out += ((tx,txc.get));
      }
      //}
    }
    reportln("TX Map Done.","progress")

    return out;
  }


  def mkSubDelimString(x : Seq[Seq[String]], delims : Seq[String]) : String = {
    x.map(_.mkString(delims(1))).mkString(delims(0));
  }
  def mkSubDelimList(x : Seq[Seq[String]], delims : Seq[String]) : java.util.List[String] = {
    x.map(_.mkString(delims(1))).toList;
  }
  
  def getMutString(ref : String, alt : String, pos : Int, 
                      getPosString : (Int => String) = {(a) => a.toString},
                      swapStrand : Boolean = false ) : String = {
    def getOrd(a : String, b : String) : (String,String) = if(swapStrand) (b,a) else (a,b);
    def getSeq(s : String) : (String) = if(swapStrand) complementSequenceString(s.reverse).mkString("") else s;
    
    if(ref.length == 1 && alt.length == 1){
      return ""+getPosString(pos)+getSeq(ref)+">"+getSeq(alt);
    } else if(ref.length == 2 && alt.length == 1 && ref.substring(0,1) == alt){
      //properly formatted single-base DEL:
      return ""+getPosString(pos+1)+"del"
    } else if(ref.length > 2 && alt.length == 1 && ref.substring(0,1) == alt){
      //properly formatted multibase DEL:
      val (a,b) = getOrd(getPosString(pos+1),getPosString(pos+ref.length-1));
      return ""+a+"_"+b+"del"
    } else if(ref.length == 1 && alt.length > 1 && alt.substring(0,1) == ref){
      //properly formatted INS:
      val (a,b) = getOrd(getPosString(pos),getPosString(pos+1));
      return ""+a+"_"+b+"ins"+getSeq(alt.substring(1));
    } else if(ref.length == 1 && alt.length > 1 && alt.substring(0,1) != ref){
      //return "g."+pos+"_"+pos+"delins"+alt;
      return ""+getPosString(pos)+"delins"+getSeq(alt);
    } else if(ref.length > 1 && alt.length == 1 && ref.substring(0,1) != alt){
      val (a,b) = getOrd(getPosString(pos),getPosString(pos+ref.length-1));
      return ""+a+"_"+b+"delins"+getSeq(alt);
    } else if(ref.length > 1 && alt.length > 1){
      //Fix oddball formatting:
      if(ref.head == alt.head){
        return getMutString(ref.tail, alt.tail, pos + 1, getPosString, swapStrand);
      } else if(ref.last == alt.last){
        return getMutString(ref.init, alt.init, pos, getPosString, swapStrand);
      } else {
        val (a,b) = getOrd(getPosString(pos),getPosString(pos+ref.length-1));
        return ""+getPosString(pos)+"_"+getPosString(pos+ref.length-1)+"delins"+getSeq(alt)
      }
    }
    //IMPOSSIBLE STATE? Malformed VCF?
    return "???";
  }
  
  def getShiftedMutString(ref : String, alt : String, pos : Int,
                      getPosString : (Int => String) = {(a) => a.toString},
                      swapStrand : Boolean = false,
                      getFlankingSeq : (Int => Option[Char])) : String = {
    def getOrd(a : String, b : String) : (String,String) = if(swapStrand) (b,a) else (a,b);
    def getSeq(s : String) : (String) = if(swapStrand) complementSequenceString(s.reverse).mkString("") else s;
    
    val refS = getSeq(ref);
    val altS = getSeq(alt);
    
    if(ref == alt){
      return "="
    } else if(ref.length > 0 && alt.length > 0 && ref.head == alt.head){
      return getShiftedMutString(ref.tail,alt.tail,pos = pos+1,getPosString = getPosString, swapStrand = swapStrand, getFlankingSeq = (i => getFlankingSeq(i+1)));
    } else if(ref.length > 0 && alt.length > 0 && ref.last == alt.last){
      return getShiftedMutString(ref.init,alt.init,pos = pos,getPosString = getPosString, swapStrand = swapStrand, getFlankingSeq = getFlankingSeq);
    } else if(ref.length == 1 && alt.length == 1){
      return ""+getPosString(pos)+getSeq(ref)+">"+getSeq(alt);
    } else if(ref.length == 1 && alt.length > 1){
        //return "g."+pos+"_"+pos+"delins"+alt;
        return ""+getPosString(pos)+"delins"+altS;
    } else if(ref.length > 1 && alt.length > 1){
        val (a,b) = getOrd(getPosString(pos),getPosString(pos+ref.length-1));
        return ""+a+"_"+b+"delins"+altS
    } else {
      //insertion or deletion: attempt realign to satisfy 3' rule:
      
      val variantLen = refS.length;
      val refSeq  = refS.toStream.map(Some(_)) ++ Stream.from(refS.length).map(getFlankingSeq(_));
      val altSeq  = altS.toStream.map(Some(_)) ++ Stream.from(refS.length).map(getFlankingSeq(_));
      val offset = refSeq.zip(altSeq).indexWhere{case (r,a) => a.isEmpty || r.isEmpty || r.get != a.get}
      
      val fref = refSeq.slice(offset,offset+refS.length).map(_.get).mkString("");
      val falt = altSeq.slice(offset,offset+altS.length).map(_.get).mkString("");
      //val varLen = math.max(finalRef.length,finalAlt.length);
      
      val offsetGS = if(! swapStrand) offset else -offset;
      
      if(ref.length == 1 && alt.length == 0){
        //single-base DEL:
        return getPosString(pos+offset)+"del"+fref
      } else if(ref.length > 1 && alt.length == 0){
        //multibase DEL:
        val (a,b) = getOrd(getPosString(pos + offset),getPosString(pos + ref.length + offset - 1));
        return ""+a+"_"+b+"del"+fref;
      } else if(ref.length == 0 && alt.length > 0){
        //insertion:
        val (a,b) = getOrd(getPosString(pos + offset - 1),getPosString(pos + offset));
        return ""+a+"_"+b+"ins"+falt
      }
    }
    /* OLD, overly complex version:
    } else if(ref.length == 1 && alt.length == 0){
      //properly formatted single-base DEL:
      val posOffset = Iterator.from(1).indexWhere{x => {
        val flank = getFlankingSeq(x)
        (! flank.isDefined) || (flank.get != refS.head)
      }}
      return ""+getPosString(pos+posOffset)+"del"
      
      //val variantLen = math.max(refS.length,altS.length);
      //val refChars  = Range(0,refS.length + variantLen).map(getFlankingSeq(_));
      //val altChars  = altS.toVector.map(Some(_)) ++ Range(altS.length,refS.length + variantLen).map(getFlankingSeq(_));
      //val offset = Iterator.from(0).map(getFlankingSeq(_)).sliding(refS.length + variantLen,1)
      
    } else if(ref.length > 1 && alt.length == 0){
      val variantLen = refS.length;
      val refSeq  = Stream.from(0).map(getFlankingSeq(_));
      val altSeq  = Stream.from(refS.length).map(getFlankingSeq(_));
      val offsetAlts = Stream.from(0).map{i => {
        (refSeq.take(i) ++ refSeq.slice(i+refS.length,i+refS.length+variantLen)).toVector
      }}
      val offsetMatch = Stream.from(0).map(i => altSeq.slice(0,i + variantLen).toVector).zip(offsetAlts).zipWithIndex.map{case ((base,offset),i) => {
        val (offStr,baseStr) = (offset.map(_.getOrElse('?')).mkString(""),base.map(_.getOrElse('?')).mkString(""))
        (offStr,baseStr,offStr(i) != '?' && offStr==baseStr,i)
      }}
      val offsetCalc = offsetMatch.grouped(variantLen).takeWhile(_.exists{case (off,base,bool,idx) => {
        bool
      }}).toVector.flatten.filter{case (off,base,bool,idx) => {
        bool
      }}
      val offset = offsetCalc.lastOption match {
        case Some(x) => x._4;
        case None => 0;
      }
      
      val (a,b) = if(swapStrand){
        (getPosString(pos+ref.length-offset-1),getPosString(pos - offset))
      } else {
        (getPosString(pos + offset),getPosString(pos+ref.length+offset-1))
      }
      return ""+a+"_"+b+"del"
      
      
      //simpler:

      
      
      //println( refSeq.slice(0,offset+variantLen+refS.length).map(_.get).mkString("") )
      //println( (repString("*",refS.length).toVector ++ refSeq.slice(refS.length,offset+variantLen+refS.length).map(_.get)).mkString("") )
      //println( (refSeq.slice(0,offset).map(_.get) ++ repString("*",refS.length).toVector ++ refSeq.slice(offset+refS.length,offset+variantLen+refS.length).map(_.get)).mkString("") )
      //println( (refSeq.slice(0,offset) ++ refSeq.slice(offset+refS.length,offset+variantLen+refS.length)).map(_.get).mkString("") )
      //println( refSeq.slice(refS.length,offset+variantLen+refS.length).map(_.get).mkString("") )
    
    } else if(ref.length == 0 && alt.length > 0){
      //properly formatted INS:
      val variantLen = math.max(refS.length,altS.length);
      val refSeq  = Stream.from(0).map(getFlankingSeq(_));
      val altSeq  = altS.map(Some(_)).toStream ++ Stream.from(refS.length).map(getFlankingSeq(_));
      val altStrings = Stream.from(0).map(i => altSeq.slice(0,i + variantLen + altS.length).toVector)
      val altInsertions = altStrings.zipWithIndex.map{case (a,i) => {
        a.slice(i,i+altS.length);
      }}
      val offsetAlts = Stream.from(0).map{i => {
        (refSeq.take(i) ++ altInsertions(i) ++ refSeq.slice(i+refS.length,i+refS.length+variantLen)).toVector
      }}
      val offsetMatch = altStrings.zip(offsetAlts).zipWithIndex.map{case ((base,offset),i) => {
        val (offStr,baseStr) = (offset.map(_.getOrElse('?')).mkString(""),base.map(_.getOrElse('?')).mkString(""))
        (offStr,baseStr,(! altInsertions(i).contains('?')) && offStr(i) != '?' && offStr==baseStr,i)
      }}
      val offsetCalc = offsetMatch.grouped(variantLen).takeWhile(_.exists{case (off,base,bool,idx) => {
        bool
      }}).toVector.flatten.filter{case (off,base,bool,idx) => {
        bool
      }}
      val offset = offsetCalc.lastOption match {
        case Some(x) => x._4;
        case None => 0;
      }
      
      val finalAlt = altInsertions(offset).map(_.get).mkString("");
      val (a,b) = if(swapStrand){
        (getPosString(pos-offset-1),getPosString(pos - offset))
      } else {
        (getPosString(pos+offset-1),getPosString(pos+offset))
      }
      return ""+a+"_"+b+"ins"+finalAlt
      
      //val (a,b) = getOrd(getPosString(pos),getPosString(pos+1));
      //return ""+a+"_"+b+"ins"+getSeq(alt.substring(1));
    } else if(ref.length == 1 && alt.length > 1 && alt.substring(0,1) != ref){
      //return "g."+pos+"_"+pos+"delins"+alt;
      return ""+getPosString(pos)+"delins"+getSeq(alt);
    } else if(ref.length > 1 && alt.length == 1 && ref.substring(0,1) != alt){
      val (a,b) = getOrd(getPosString(pos),getPosString(pos+ref.length-1));
      return ""+a+"_"+b+"delins"+getSeq(alt);
    } else if(ref.length > 1 && alt.length > 1){
      //Fix oddball formatting:
      val (a,b) = getOrd(getPosString(pos),getPosString(pos+ref.length-1));
      return ""+getPosString(pos)+"_"+getPosString(pos+ref.length-1)+"delins"+getSeq(alt)
    }*/
    //IMPOSSIBLE STATE? Malformed VCF?
    return "???";
  }
  

  def buildTXUtilFromString(s : String) : TXUtil = {
    val cells = s.split("\t"); 
    
    return new TXUtil(
           geneID = cells(0), 
           txID = cells(1), 
           chrom = cells(2), 
           strand = cells(3).charAt(0), 
           gSpans = cells(4).split(";").map(pair => {
             val pairSplit = pair.split(",");
             (string2int(pairSplit(0)), string2int(pairSplit(1)));
           }).toVector, 
           gCdsStart = string2int(cells(5)), 
           gCdsEnd = string2int(cells(6)),
           seqGS = cells(7).toVector,
           leadingExonSeq = cells(8).split(";").toVector,
           trailingExonSeq = cells(9).split(";").toVector,
           gCdsNS = (string2int(cells(10).split(";")(0)), string2int(cells(10).split(";")(1)))
    )
  }
  
}