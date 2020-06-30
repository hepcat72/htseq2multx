#!/bin/bash

#USAGE: ./run_tests.sh [0|1|2] [last_test_number]
#First argument is a debug value. 1 adds --debug to htseq2multx.pl commands to
#cause the script to echo the fastq-multx command.  2 also causes the script to
#stop upon the end of the first failed test.  3 echoes all the test script's
#commands

TEST_DATA=test_data
TEST_OUTPUT=test_output
BARCODE_SPLITTER=../htseq2multx.pl
DEBUG=0
DEBUG_OPT=""
if [ $# -ne 0 ] && [ $1 -ne 0 ]; then
  DEBUG=$1
  DEBUG_OPT="--debug"
fi
LAST_TEST=""
if [ $# -gt 1 ]; then
  LAST_TEST="test_$2"
fi
rm -rf "${TEST_OUTPUT}"
mkdir -p "${TEST_OUTPUT}"

#Perl code to allow fastq-multx debug message to not cause a "difference"
PERL_DEBUG_FILTER='$c=0;$d="";while(<STDIN>){if($_ !~ /#.*fastq-multx|^1d0$/){print;$c++}elsif(/#.*fastq-multx/){$d.=$_}}if($c && $ARGV[0]){$d =~ s/< #/#/g;print "\033[1;94m$d\033[0m"}exit($c)'

EXIT_STATUS=0

function check_test_results() {
    if [ $EXIT_CODE -ne $EXPECTED_EXIT_CODE ];then
        EXIT_STATUS=1
        echo -e "\033[1;31mFAILED ${TEST} EXIT($EXIT_CODE) != EXPECTED($EXPECTED_EXIT_CODE)\033[0m"
    else
        echo -e "\033[1;32mPASSED\033[0m ${TEST} EXIT($EXIT_CODE)"
    fi
    for OUTPUT in $OUTPUTS; do
        filename=$(basename "$OUTPUT")
        extension="${filename##*.}"
        op=''
        if [[ $extension == 'gz' ]]; then
            diff <(gzip -dc "$TEST_OUTPUT/${TEST}_${OUTPUT}") <(gzip -dc "${TEST_DATA}/${TEST}_${OUTPUT}") | perl -e "$PERL_DEBUG_FILTER" $DEBUG
        else
            diff "${TEST_OUTPUT}/${TEST}_${OUTPUT}" "${TEST_DATA}/${TEST}_${OUTPUT}" | perl -e "$PERL_DEBUG_FILTER" $DEBUG
        fi
        if [ $? -ne 0 ]; then
            EXIT_STATUS=1
            echo -e "\033[1;31mFAILED ${TEST}_${OUTPUT}\033[0m"
        else
            echo -e "\033[1;32mPASSED\033[0m ${TEST}_${OUTPUT}"
        fi
    done
    for NOT_OUTPUT in $NOT_OUTPUTS; do
        if [ -f "$NOT_OUTPUT" ]
        then
            EXIT_STATUS=1
            echo -e "\033[1;31mFAILED ${TEST}_${NOT_OUTPUT}\033[0m"
        else
            echo -e "\033[1;32mPASSED\033[0m ${TEST}_${NOT_OUTPUT}"
        fi
    done
    if [ $DEBUG -gt 1 ] && [ $EXIT_STATUS -gt 0 ]; then
      exit 1
    fi
    if [ "$LAST_TEST" == "$TEST" ]; then
      exit 0
    fi
}

NOT_OUTPUTS=""

if [ $DEBUG -gt 2 ]; then
  set -x
fi

# Test split file output
TEST=test_1
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test separate index file
TEST=test_2
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out BC1-read-2.out BC2-read-2.out BC3-read-2.out BC4-read-2.out unmatched-read-2.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" "${TEST_DATA}/barcode_splitter_index.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 2 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test special characters in sample names (utf-8)
TEST=test_3
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC_1-read-1.out BC_2-read-1.out BC_3-read-1.out BC_4-read-1.out unmatched-read-1.out BC_1-read-2.out BC_2-read-2.out BC_3-read-2.out BC_4-read-2.out unmatched-read-2.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes_odd_sample_names.txt" "${TEST_DATA}/barcode_splitter1.fastq" "${TEST_DATA}/barcode_splitter_index.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 2 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test variable length barcodes
TEST=test_4
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes_vary_length.txt" "${TEST_DATA}/barcode_splitter1.fastq" "${TEST_DATA}/barcode_splitter_index.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 2 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test gzip input and output
TEST=test_5
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out.gz BC2-read-1.out.gz BC3-read-1.out.gz BC4-read-1.out.gz unmatched-read-1.out.gz summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq.gz" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --gzipout --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test dual barcodes
TEST=test_6
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.fastq BC1-read-2.fastq BC1-read-3.fastq BC2-read-1.fastq BC2-read-2.fastq BC2-read-3.fastq BC3-read-1.fastq BC3-read-2.fastq BC3-read-3.fastq BC4-read-1.fastq BC4-read-2.fastq BC4-read-3.fastq error.out summary.out unmatched-read-1.fastq unmatched-read-2.fastq unmatched-read-3.fastq"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --bcfile "${TEST_DATA}/barcode_splitter_barcodes_dual.txt" --prefix "${TEST_OUTPUT}/${TEST}_" "${TEST_DATA}/barcode_splitter1.fastq" "${TEST_DATA}/barcode_splitter_index.fastq" "${TEST_DATA}/barcode_splitter_index_2.fastq" --idxread 2 3 --mismatches 1 --split_all $DEBUG_OPT 1> "${TEST_OUTPUT}/${TEST}_summary.out" 2> "${TEST_OUTPUT}/${TEST}_error.out"
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test old Illumina style fastq ids
TEST=test_7
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out BC1-read-2.out BC2-read-2.out BC3-read-2.out BC4-read-2.out unmatched-read-2.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1_oldstyle_ids.fastq" "${TEST_DATA}/barcode_splitter_index_oldstyle_ids.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 2 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test new Illumina style fastq ids
TEST=test_8
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out BC1-read-2.out BC2-read-2.out BC3-read-2.out BC4-read-2.out unmatched-read-2.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1_newstyle_ids.fastq" "${TEST_DATA}/barcode_splitter_index_newstyle_ids.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 2 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test no splitting of index file by default
TEST=test_9
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out summary.out error.out"
NOT_OUTPUTS="BC1-read-2.out BC2-read-2.out BC3-read-2.out BC4-read-2.out unmatched-read-2.out"
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" "${TEST_DATA}/barcode_splitter_index.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 2 $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test split all files by default when all are index files
TEST=test_10
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test that sequences shorter than barcodes are treated as unmatched
TEST=test_11
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 0 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1_short.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test output sequence file does not exist
TEST=test_12
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "baddir/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=1
check_test_results


# Test input sequence file does not exist
TEST=test_13
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "baddir/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=10
check_test_results


# Test input barcode file does not exist
TEST=test_14
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "baddir/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=5
check_test_results


# Test output unmatched file exists
TEST=test_15
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
touch ${TEST_OUTPUT}/${TEST}_unmatched-read-1.out
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=9
check_test_results


# Test output sample file exists
TEST=test_16
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
touch ${TEST_OUTPUT}/${TEST}_BC1-read-1.out
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=9
check_test_results


# Test gzipped output sequence file does not exist
TEST=test_17
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "baddir/${TEST}_" --suffix .out --idxread 1 --split_all --gzipout $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=1
check_test_results


# Test gzipped input sequence file does not exist
TEST=test_18
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
# NOTE: Skipping check of error.out since gzip messages can differ bewteen systems
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "baddir/barcode_splitter1.fastq.gz" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=10
check_test_results


# Test gzipped output unmatched file exists
TEST=test_19
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
touch ${TEST_OUTPUT}/${TEST}_unmatched-read-1.out.gz
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all --gzipout $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=9
check_test_results


# Test gzipped output sample file exists
TEST=test_20
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS=""
touch ${TEST_OUTPUT}/${TEST}_BC1-read-1.out.gz
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/barcode_splitter1.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all --gzipout $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=9
check_test_results


# Test error for special characters in sample names (utf-8)
TEST=test_21
echo -e "\033[1;33mSKIPPING\033[0m $TEST"


# Test that version in barcode splitter and setup are the same (and that thus,
# single source versioning is working)
TEST=test_22
echo -e "\033[1;33mSKIPPING\033[0m $TEST"


# Test split file output
TEST=test_23
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out summary.out error.out"
NOT_OUTPUTS=""
mkfifo "${TEST_DATA}/${TEST}_read1_fifo.fastq"
cat "${TEST_DATA}/barcode_splitter1.fastq" > "${TEST_DATA}/${TEST}_read1_fifo.fastq" &
$BARCODE_SPLITTER --mismatches 2 --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" "${TEST_DATA}/${TEST}_read1_fifo.fastq" --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
rm "${TEST_DATA}/${TEST}_read1_fifo.fastq"
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test duplicate barcode rows
TEST=test_24
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="summary.out error.out"
NOT_OUTPUTS="BC1-read-1.out BC2-read-1.out BC3-read-1.out BC4-read-1.out unmatched-read-1.out"
$BARCODE_SPLITTER --bcfile "${TEST_DATA}/barcode_splitter_barcodes_dupes.txt" --prefix "${TEST_OUTPUT}/${TEST}_" "${TEST_DATA}/barcode_splitter1.fastq" "${TEST_DATA}/barcode_splitter_index.fastq" "${TEST_DATA}/barcode_splitter_index_2.fastq" --idxread 2 3 $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=12
check_test_results


# Test gzip input and output bash process substitution
TEST=test_25
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out.gz BC2-read-1.out.gz BC3-read-1.out.gz BC4-read-1.out.gz unmatched-read-1.out.gz summary.out error.out"
NOT_OUTPUTS=""
$BARCODE_SPLITTER --mismatches 2 --gzipin --gzipout --bcfile "${TEST_DATA}/barcode_splitter_barcodes.txt" <(cat "${TEST_DATA}/barcode_splitter1.fastq.gz") --prefix "${TEST_OUTPUT}/${TEST}_" --suffix .out --idxread 1 --split_all $DEBUG_OPT 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


# Test the fastq-multx executable itself and its handling of gzip data
TEST=test_26
rm -f ${TEST_OUTPUT}/${TEST}_* 2> /dev/null
OUTPUTS="BC1-read-1.out.gz BC2-read-1.out.gz BC3-read-1.out.gz BC4-read-1.out.gz unmatched-read-1.out.gz summary.out error.out"
NOT_OUTPUTS=""
fastq-multx -x -d 1 -B "${TEST_DATA}/barcode_splitter_barcodes_fqmx.txt" -m 2 -M 2 -b "${TEST_DATA}/barcode_splitter1.fastq.gz" -o "${TEST_OUTPUT}/${TEST}_%-read-1.out.gz" 2> ${TEST_OUTPUT}/${TEST}_error.out 1> ${TEST_OUTPUT}/${TEST}_summary.out
EXIT_CODE=$?
EXPECTED_EXIT_CODE=0
check_test_results


exit $EXIT_STATUS
