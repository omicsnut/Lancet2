#include "lancet/hts/bgzf_ostream.h"

#include "lancet/base/logging.h"

extern "C" {
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
}

#include <filesystem>
#include <ios>

#include <cstdio>
#include <cstdlib>

namespace lancet::hts {

namespace detail {

auto BgzfStreambuf::Open(std::filesystem::path const& path, char const* mode) -> bool {
  if (mFilePtr != nullptr) Close();

  mFileName = path;
  mFilePtr = bgzf_open(mFileName.c_str(), mode);
  return mFilePtr != nullptr;
}

void BgzfStreambuf::Close() {
  if (mFilePtr != nullptr) {
    // bgzf_close on a network stream (S3/GCS) blocks until all buffered multipart fragments
    // are uploaded. If it returns -1 (authentication failure, network error) and we ignore it,
    // the pipeline terminates without a diagnostic message.
    if (bgzf_close(mFilePtr) < 0) {
      LOG_CRITICAL("Failed to close BGZF stream. If writing to cloud URIs, check remote "
                   "auth/permissions: {}",
                   mFileName.string())
      std::exit(EXIT_FAILURE);
    }
    mFilePtr = nullptr;
  }
}

auto BgzfStreambuf::uflow() -> int {
  if (mCurrPos != SENTINEL_BUFFER_POSITION) {
    auto const res = mCurrPos;
    mCurrPos = SENTINEL_BUFFER_POSITION;
    return res;
  }

  if (mFilePtr == nullptr) return EOF;

  mCurrPos = bgzf_getc(mFilePtr);
  switch (mCurrPos) {
    case -1:
    case -2:
      return EOF;
    default:
      return mCurrPos;
  }
}

auto BgzfStreambuf::underflow() -> int {
  if (mFilePtr == nullptr) return EOF;
  if (mCurrPos != SENTINEL_BUFFER_POSITION) return mCurrPos;

  mCurrPos = bgzf_getc(mFilePtr);
  switch (mCurrPos) {
    case -1:
    case -2:
      return EOF;
    default:
      return mCurrPos;
  }
}

auto BgzfStreambuf::overflow(int dat) -> int {
  auto const cdat = static_cast<char>(dat);
  auto const num_bytes = bgzf_write(mFilePtr, &cdat, 1);
  return num_bytes < 0 ? static_cast<int>(num_bytes) : dat;
}

auto BgzfStreambuf::xsputn(char const* data, std::streamsize len) -> std::streamsize {
  if (mFilePtr == nullptr) return 0;
  return bgzf_write(mFilePtr, data, static_cast<std::size_t>(len));
}

auto BgzfStreambuf::sync() -> int {
  if (mFilePtr == nullptr) {
    return 0;
  }

  // std::ostream::flush invokes sync(). HTSlib manages its own buffered remote streams,
  // so we must forward the flush to bgzf_flush to ensure data reaches the network layer.
  if (bgzf_flush(mFilePtr) < 0) {
    LOG_CRITICAL(
        "Failed to flush BGZF stream. If writing to cloud URIs, check remote auth/permissions: {}",
        mFileName.string())
    std::exit(EXIT_FAILURE);
  }
  return 0;
}

}  // namespace detail

auto BgzfOstream::Open(std::filesystem::path const& path, BgzfFormat ofmt) -> bool {
  mOutFmt = ofmt;
  auto result = mBgzfBuffer.Open(path, "w");
  rdbuf(&mBgzfBuffer);
  return result;
}

void BgzfOstream::Close() {
  mBgzfBuffer.Close();
  if (mOutFmt != BgzfFormat::UNSPECIFIED) BuildIndex();
}

void BgzfOstream::BuildIndex() {
  int result = 0;
  switch (mOutFmt) {
    case BgzfFormat::VCF:
      result = tbx_index_build(mBgzfBuffer.mFileName.c_str(), 0, &tbx_conf_vcf);
      break;
    case BgzfFormat::GFF:
      result = tbx_index_build(mBgzfBuffer.mFileName.c_str(), 0, &tbx_conf_gff);
      break;
    case BgzfFormat::BED:
      result = tbx_index_build(mBgzfBuffer.mFileName.c_str(), 0, &tbx_conf_bed);
      break;
    default:
      break;
  }

  if (result < 0) {
    LOG_CRITICAL("Failed to build secondary tabix index for output file: {}",
                 mBgzfBuffer.mFileName.string())
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace lancet::hts
