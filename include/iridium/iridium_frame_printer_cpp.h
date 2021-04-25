/* -*- c++ -*- */
/*
 * Copyright 2021 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_IRIDIUM_IRIDIUM_FRAME_PRINTER_CPP_H
#define INCLUDED_IRIDIUM_IRIDIUM_FRAME_PRINTER_CPP_H

#include <gnuradio/sync_block.h>
#include <iridium/api.h>

namespace gr {
namespace iridium {

/*!
 * \brief <+description of block+>
 * \ingroup iridium
 *
 */
class IRIDIUM_API iridium_frame_printer_cpp : virtual public gr::sync_block
{
public:
    typedef std::shared_ptr<iridium_frame_printer_cpp> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of iridium::iridium_frame_printer_cpp.
     *
     * To avoid accidental use of raw pointers, iridium::iridium_frame_printer_cpp's
     * constructor is in a private implementation
     * class. iridium::iridium_frame_printer_cpp::make is the public interface for
     * creating new instances.
     */
    static sptr make(std::string file_info);
};

} // namespace iridium
} // namespace gr

#endif /* INCLUDED_IRIDIUM_IRIDIUM_FRAME_PRINTER_CPP_H */
