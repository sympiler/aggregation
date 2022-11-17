/*
 * exceptions.h
 *
 *  Created on: Sep. 27, 2020
 *      Author: shujianqian
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_


#include <exception>
#include <string>

 namespace sym_lib {

 class missing_arg_error : public std::runtime_error
 {
 public:
  missing_arg_error (std::string arg, std::string msg="Argument missing")
   : std::runtime_error(msg)
   {
    arg_ = arg;
   }

  std::string arg() const { return arg_; }

 private:
  std::string arg_;
 };

 class open_file_error : public std::runtime_error
 {
 public:
  open_file_error (std::string filename, std::string msg="Failed to open file")
  : std::runtime_error(msg)
  {
   filename_ = filename;
  }

  std::string filename() const { return filename_; }

 private:
  std::string filename_;
 };

 class read_file_error : public open_file_error
  {
  public:
   read_file_error (std::string filename, std::string msg="Failed to read file")
   : open_file_error(filename, msg)
   {}
  };

 class write_file_error : public open_file_error
 {
 public:
  write_file_error (
    std::string filename,
    std::string msg="Failed to write to file"
  )
  : open_file_error(filename, msg)
  {}
 };

 class mtx_error : public std::runtime_error
 {
 public:
  mtx_error (std::string filename, std::string msg="Error loading matrix")
  : std::runtime_error(msg)
  {
   filename_ = filename;
  }

  std::string filename() const { return filename_; }

 private:
  std::string filename_;
 };

 class mtx_header_error : public mtx_error
 {
 public:
  mtx_header_error (
    std::string filename="Unknown",
    std::string msg="Invalid matrix header"
  )
  : mtx_error(filename, msg)
  {}
 };

 class mtx_format_error : public mtx_error
 {
 public:
  mtx_format_error
  (
    std::string expected_format,
    std::string got_format,
    std::string filename = "Unknown",
    std::string msg = "Matrix format mismatch"
  )
  : mtx_error(filename, msg)
  {
   expected_format_ = expected_format;
   got_format_ = got_format;
  }

  std::string expected_format() const { return expected_format_; }
  std::string got_format() const { return got_format_; }

 private:
  std::string expected_format_;
  std::string got_format_;
 };

 class mtx_arith_error : public mtx_error
 {
 public:
  mtx_arith_error
  (
    std::string expected_arith,
    std::string got_arith,
    std::string filename="Unknown",
    std::string msg="Matrix arithmetic mismatch"
  )
  : mtx_error(filename, msg)
  {
   expected_arith_ = expected_arith;
   got_arith_ = got_arith;
  }
  std::string expected_arith() const { return expected_arith_; }
  std::string got_arith() const { return got_arith_; }

 private:
  std::string expected_arith_;
  std::string got_arith_;
 };

}  // namespace format

#endif /* EXCEPTIONS_H_ */
