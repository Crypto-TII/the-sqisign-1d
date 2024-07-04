// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rng.h>
#include <sig.h>
#include <api.h>

// Some hardcoded vectors for profiling

#if CRYPTO_PUBLICKEYBYTES == 64
const unsigned char pk_lvl1[] = {0x77,0x69,0xC5,0x34,0xCF,0xA4,0x02,0x68,0x53,0xF0,0x41,0xF1,0x4A,0xB3,0xAA,0x97,0x5E,0x2D,0x15,0xC4,0xE5,0x14,0x73,0xC5,0xEC,0x44,0x81,0x50,0x9A,0x09,0x58,0x1C,0x2A,0x27,0x46,0x9E,0x49,0x10,0x69,0xC9,0xB3,0x60,0x9B,0x4E,0x3B,0xF3,0x1C,0xFB,0xF2,0x26,0xBF,0xF0,0xA4,0x80,0xBC,0x48,0x43,0x80,0xD8,0x22,0x00,0x98,0x13,0x02}; 
const unsigned char sk_lvl1[] = {0x77,0x69,0xC5,0x34,0xCF,0xA4,0x02,0x68,0x53,0xF0,0x41,0xF1,0x4A,0xB3,0xAA,0x97,0x5E,0x2D,0x15,0xC4,0xE5,0x14,0x73,0xC5,0xEC,0x44,0x81,0x50,0x9A,0x09,0x58,0x1C,0x2A,0x27,0x46,0x9E,0x49,0x10,0x69,0xC9,0xB3,0x60,0x9B,0x4E,0x3B,0xF3,0x1C,0xFB,0xF2,0x26,0xBF,0xF0,0xA4,0x80,0xBC,0x48,0x43,0x80,0xD8,0x22,0x00,0x98,0x13,0x02,0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xAF,0x71,0x87,0xCD,0xA3,0x13,0xF3,0x7F,0xF9,0xF3,0xFD,0xAB,0xC6,0x25,0x48,0xF9,0x69,0xF2,0x5C,0xAB,0x3B,0xF9,0xAD,0x2C,0x79,0x07,0x22,0x7C,0xA6,0xAD,0x10,0xF3,0x64,0x9A,0xD7,0xF3,0x68,0xB6,0xCD,0x97,0xD7,0x6F,0xF5,0x67,0x29,0xE9,0x03,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xFA,0xB3,0x43,0xC8,0x0E,0xF4,0x16,0xB8,0xBE,0x87,0x4E,0x2F,0x6D,0x2D,0xA7,0x38,0x00,0x81,0x24,0x17,0x38,0x44,0xDB,0x60,0x19,0xB7,0xB8,0x23,0xC4,0x9E,0xB7,0xD1,0x1B,0x46,0x6C,0x99,0xDC,0xCD,0xF9,0xBD,0x34,0x1F,0x73,0x0E,0x92,0x32,0xFE,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xEE,0x8D,0x38,0xC1,0x1F,0xF3,0xB5,0x0E,0x8B,0xED,0x00,0x09,0x48,0xC0,0x60,0x4D,0x63,0xCB,0xC8,0x94,0x3C,0xAD,0x5A,0x0F,0xCD,0x1A,0xBF,0x81,0x30,0x81,0xFB,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xBF,0xC9,0x2E,0x58,0x2C,0xC3,0x50,0x79,0xFA,0x93,0x64,0x3B,0x1F,0x12,0xB3,0xF0,0x86,0x43,0xF6,0xBB,0xAA,0x1B,0xB4,0x3E,0xBB,0x60,0x1D,0x25,0x4D,0x03,0x0C,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x22,0x44,0xAD,0x1D,0x60,0x35,0xD2,0x53,0x02,0x92,0xE1,0xB9,0xEB,0x61,0x10,0xAE,0x11,0xD9,0xB2,0x0A,0xD3,0xFC,0xF0,0x50,0x49,0xCC,0x5E,0x1B,0x8F,0x77,0x82,0x0C,0xA9,0xEB,0x43,0x9F,0x20,0xFF,0x06,0x17,0x2E,0x1C,0x13,0xDE,0xDB,0x33,0xE7,0x24,0xD2,0x60,0x90,0xA5,0x69,0xAB,0xE3,0x5E,0xE8,0x37,0xEE,0x18,0x23,0x15,0x0B,0x15,0xB6,0xF1,0x65,0xB3,0x24,0x08,0x65,0x24,0xB7,0x41,0xFA,0x46,0x88,0x38,0x04,0x63,0x88,0x94,0x92,0x10,0xE1,0x4C,0xCA,0x16,0x7B,0x0A,0x02,0x5D,0x61,0xDE,0xBE,0x33,0x8B,0x96,0x68,0x43,0xBF,0xBA,0xE5,0x4D,0x78,0x96,0x3E,0x8D,0xE6,0x27,0x34,0x92,0xC5,0x32,0x7E,0x2A,0x29,0x4A,0xA4,0x63,0xD0,0x90,0x05,0xE6,0x33,0x61,0x2C,0x1E,0x45,0xB1,0x63,0x5A,0xFA,0x8F,0xEB,0x9D,0x5F,0x44,0x6B,0x43,0xFB,0x62,0x1A,0x03,0x91,0xE8,0x86,0x7B,0xF3,0x95,0x4A,0xE6,0x16,0xAC,0x53,0x64,0xED,0x38,0x74,0x08,0x78,0x6F,0xC0,0x3C,0x45,0x5B,0xEE,0x6C,0x80,0x9C,0xFF,0xB6,0x5D,0x35,0x20,0xF6,0xDC,0x2E,0x27,0x75,0xAE,0x34,0x72,0x5D,0x32,0xEA,0x49,0x0F,0x9D,0xFE,0x0D,0x30,0x28,0x74,0xF2,0xB1,0x85,0xB1,0xF7,0xD7,0x83,0x20,0x20,0xBA,0x97,0xD6,0x50,0xBB,0x09,0x22,0x17,0xF4,0x57,0x6E,0xD8,0xC3,0xEF,0x8D,0x95,0x85,0xD4,0x0A,0x56,0x13,0xF5,0xD1,0x32,0x31,0xEE,0x17,0x56,0x6C,0x49,0x97,0xCD,0x47,0xFE,0x11,0x3B,0x95,0x6B,0xC9,0x54,0x0B,0x68,0x5B,0xB1,0x03,0xD6,0xDD,0xCA,0x75,0x70,0x80,0xBE,0x19,0x00,0xD7,0x1D,0x44,0xC0,0xCE,0x81,0x0E,0x0D,0xF1,0xA0,0x17,0x25,0xF4,0x70,0xE7,0x7B,0xE8,0xAA,0xB8,0x14,0xDB,0x3E,0x51,0xB8,0x6D,0xCF,0xF5,0xAD,0xD8,0x28,0x1B,0xD9,0xF9,0xBC,0x67,0x04,0xB8,0x46,0xE3,0x38,0xEB,0xFD,0xE4,0x3B,0x72,0x44,0x28,0x05,0x46,0x14,0x73,0xB6,0xD2,0x76,0x6C,0xF9,0x30,0x74,0xE4,0x78,0x0D,0xE5,0x0B,0x36,0x1F,0x43,0xD8,0x2C,0xF2,0x83,0x8D,0x0C,0x4A,0x2A,0x00,0xFE,0x44,0xBB,0xA2,0xC1,0x06,0xC5,0x7B,0x99,0x33,0xD9,0x34,0xE4,0x56,0x23,0x88,0xAD,0xA3,0xD5,0x1C,0x90,0x3C,0x9C,0x77,0x3A,0xFB,0x3A,0x0E,0xB8,0x62,0xCF,0x41,0x0B,0xD2,0xD3,0x79,0x3D,0x4C,0x53,0x18,0xC9,0x35,0x60,0x26,0xCB,0x35,0xB4,0xCF,0x8E,0x37,0x01,0x34,0x2F,0xDA,0xA1,0x19,0xEA,0x49,0x38,0xB2,0xCC,0x71,0x5C,0xF0,0x7A,0x4C,0x9C,0xE1,0x68,0x18,0xCF,0x15,0xD3,0x71,0xED,0x6E,0x45,0x64,0xA9,0x31,0x2F,0xEB,0x57,0x0C,0xAB,0x2D,0x28,0x71,0x45,0xCF,0xAE,0x87,0xF1,0xD1,0xD5,0x6D,0x4D,0xBA,0x77,0x6F,0xA4,0x48,0xA7,0xD1,0x38,0xFE,0x3C,0xE2,0x80,0x13,0x68,0x6F,0xA0,0x5D,0xA1,0x08};
const unsigned char sm_lvl1[] = {0xAC,0xCF,0xF0,0x35,0x98,0xEF,0xB1,0x3D,0x03,0x07,0xBC,0xDF,0x97,0x9B,0x73,0x48,0xA9,0x51,0xAB,0x00,0xCF,0x92,0xF0,0x47,0x08,0xD0,0x06,0x7E,0x72,0x06,0x13,0x97,0xA8,0x1D,0x92,0xE6,0x11,0xFF,0x68,0x01,0xE1,0x4E,0x55,0xE5,0xFF,0xA9,0xBE,0xFE,0x78,0x07,0xB3,0x8A,0x9D,0x40,0x41,0x64,0xE0,0x65,0x7B,0x02,0x56,0x4E,0x2B,0x26,0x38,0x3D,0xB0,0x05,0x01,0x01,0x06,0x13,0x58,0x57,0xCE,0xFC,0xE6,0x1E,0xFD,0x00,0x9D,0x8A,0x60,0x00,0x3E,0x09,0x8F,0x43,0xDB,0x06,0x17,0x3F,0x2A,0x68,0x43,0x0A,0xC1,0xEC,0x70,0x03,0x19,0xC7,0xF7,0xB4,0x7F,0x81,0x01,0x8A,0xFB,0x00,0xA7,0x06,0x7A,0x1B,0x11,0x5D,0xBD,0x3B,0xE8,0x05,0xFD,0x2D,0x24,0x03,0x0C,0x91,0x50,0xDD,0x48,0x06,0xC4,0x8D,0xEA,0x84,0xEC,0xAC,0x2C,0x59,0x48,0x00,0x00,0x43,0x23,0x1B,0xD6,0xB0,0x8C,0x2B,0x7C,0x39,0xD9,0x32,0x06,0xE3,0xCD,0x79,0x73,0x10,0x01,0xE6,0xD5,0xBF,0x78,0xDB,0x04,0x03,0xDB,0xD2,0x05,0x4F,0x16,0x3B,0x51,0xA4,0xCF,0xB6,0x01,0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
const unsigned char msg_lvl1[] = {0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
const size_t mlen_lvl1 = 33;
#elif CRYPTO_PUBLICKEYBYTES == 96
const unsigned char pk_lvl1[] = {0xFB,0xDB,0x06,0x7F,0xED,0x9C,0xDF,0x27,0xB1,0xDC,0x94,0x4A,0xE1,0x55,0xD4,0x81,0x7F,0xE2,0x25,0x85,0xBE,0x32,0xB3,0x54,0xE6,0x2D,0x89,0x0C,0xA2,0xFF,0x79,0x39,0xD6,0x3C,0x59,0x14,0x4A,0xD6,0x6A,0x77,0x42,0x7D,0x58,0xC7,0x53,0xAD,0x2A,0x00,0xB0,0xFB,0xBA,0x29,0xCA,0x34,0x47,0x51,0x40,0xCB,0x3F,0xFF,0xFC,0x12,0x2D,0xDF,0x1F,0xEE,0x1A,0xDD,0xFE,0xD3,0x89,0x31,0x8A,0x61,0xA0,0xE9,0x98,0x14,0x7A,0x30,0x81,0xCA,0x1B,0x71,0xD1,0x3B,0xFC,0x80,0x6F,0x13,0xF7,0x02,0x07,0x19,0x2A,0x01};
const unsigned char sk_lvl1[] = {0xFB,0xDB,0x06,0x7F,0xED,0x9C,0xDF,0x27,0xB1,0xDC,0x94,0x4A,0xE1,0x55,0xD4,0x81,0x7F,0xE2,0x25,0x85,0xBE,0x32,0xB3,0x54,0xE6,0x2D,0x89,0x0C,0xA2,0xFF,0x79,0x39,0xD6,0x3C,0x59,0x14,0x4A,0xD6,0x6A,0x77,0x42,0x7D,0x58,0xC7,0x53,0xAD,0x2A,0x00,0xB0,0xFB,0xBA,0x29,0xCA,0x34,0x47,0x51,0x40,0xCB,0x3F,0xFF,0xFC,0x12,0x2D,0xDF,0x1F,0xEE,0x1A,0xDD,0xFE,0xD3,0x89,0x31,0x8A,0x61,0xA0,0xE9,0x98,0x14,0x7A,0x30,0x81,0xCA,0x1B,0x71,0xD1,0x3B,0xFC,0x80,0x6F,0x13,0xF7,0x02,0x07,0x19,0x2A,0x01,0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xA6,0x25,0xA3,0xE9,0x72,0x58,0x8F,0xA8,0x1D,0xF7,0x63,0x37,0x83,0xA0,0x7C,0x17,0x06,0x82,0x14,0x32,0xC9,0x75,0xD2,0x7F,0x15,0xA3,0x4F,0xA6,0x7B,0x9D,0x22,0x5A,0xDA,0xCE,0xD4,0xDE,0x4B,0x43,0xD8,0xFC,0x47,0xCE,0x7D,0x27,0xED,0x46,0x91,0x5D,0x23,0xB7,0xC4,0x24,0x09,0x7C,0x49,0xDA,0xF3,0x08,0x3D,0x63,0x80,0x19,0xA4,0x08,0xE2,0x60,0x0D,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x65,0x08,0x56,0xD9,0x6E,0xDE,0x15,0xB6,0x89,0xCE,0x0A,0xAC,0x09,0xEB,0xD7,0x57,0x79,0xA4,0xC3,0xE7,0x81,0x6D,0x66,0x39,0x2E,0xF0,0x0A,0xCA,0x5A,0x17,0x4C,0x39,0x1A,0xBB,0x6D,0xF1,0x70,0x6D,0x1A,0xCE,0x4E,0xC4,0x56,0xC7,0x17,0x03,0xE7,0x33,0x45,0x8A,0x34,0xE4,0xA6,0x5E,0xE8,0xBD,0x4A,0x59,0x0E,0x8F,0xB7,0x6B,0xDC,0x3F,0xA5,0x1A,0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x45,0x1F,0x15,0xE0,0x2A,0xF5,0xA3,0x79,0x65,0x77,0xD5,0x10,0x00,0xD0,0xAE,0xDA,0xAD,0x34,0x4F,0xB7,0xBD,0x70,0xCE,0x10,0x69,0x0A,0x93,0xBF,0xE3,0xEB,0xFE,0x10,0x66,0xB1,0x7B,0x6B,0x9F,0x8B,0x19,0x76,0xEC,0x07,0x4A,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xD6,0x0B,0x88,0x66,0xD7,0xE3,0x4F,0x7D,0xC0,0x7F,0xD5,0xB3,0x84,0xEA,0x6C,0xDA,0xD7,0xE5,0xC3,0x3C,0xAC,0x15,0xB2,0x03,0x5B,0xB1,0x1F,0x6B,0x63,0x7C,0x4B,0x61,0xC0,0xC1,0x9B,0xF3,0xC2,0xBB,0xD7,0x2A,0x14,0x1F,0xDC,0x01,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x33,0x8A,0xFD,0xE2,0x24,0x58,0x65,0xC8,0x94,0x07,0x19,0xA4,0x20,0x84,0xC0,0xD6,0xC6,0xEC,0x81,0x9B,0x6E,0x54,0x21,0x2F,0xE9,0x11,0x25,0x7A,0x48,0xE3,0x58,0xDB,0xAE,0x38,0x21,0x87,0xFD,0x9F,0x54,0xE9,0x8E,0x0D,0x67,0x63,0x5E,0x55,0x12,0x03,0xB6,0x8F,0xE4,0xA4,0x97,0x1C,0x9A,0x77,0x35,0x04,0x5C,0xED,0xCE,0x83,0xAD,0x8A,0xDC,0x26,0x0D,0x6F,0x65,0xFD,0x60,0x51,0xDA,0xB1,0x8F,0x00,0x61,0x1F,0xCE,0x77,0x74,0x34,0x08,0x0A,0x5B,0x36,0x67,0x61,0xF8,0x9C,0xF7,0xE2,0xD8,0x84,0xDB,0x02,0xF1,0x4C,0x4E,0xB0,0x66,0xD8,0x25,0x52,0x64,0x89,0xC8,0xD6,0xA7,0x44,0x35,0xE1,0x2B,0x63,0x32,0xAB,0x2B,0x7A,0x23,0x00,0x20,0xC1,0x9B,0xEA,0x53,0x92,0x24,0x14,0x7E,0xC5,0xB2,0x07,0xBA,0x18,0xCD,0xBA,0xCF,0xB4,0x27,0x23,0xCF,0x94,0x95,0x00,0x10,0x69,0x30,0xC9,0xB3,0x73,0x38,0x43,0x96,0xE1,0x55,0xC2,0xDB,0xA1,0x39,0x70,0xB5,0x44,0x27,0x40,0x02,0x04,0xB1,0x3B,0x6D,0x46,0x29,0xF8,0x99,0xD8,0xEC,0x29,0xA9,0x55,0x7D,0x67,0xC7,0x2E,0x1F,0xE5,0xC4,0x0A,0x04,0x65,0x2D,0x71,0x18,0x01,0xD5,0x2D,0x54,0x65,0x77,0xC4,0xFE,0xAE,0x74,0x0E,0x4C,0xDB,0x3D,0xAF,0x6D,0x62,0xF4,0xC9,0xF0,0x1B,0xBE,0xDC,0x8C,0xC2,0x81,0x4E,0xA6,0x04,0x3C,0x51,0x78,0xB1,0xB6,0x4D,0x23,0x9A,0x26,0xCC,0x5F,0xB1,0xBD,0x3B,0xBB,0xF6,0x47,0xE5,0x30,0x00,0x5F,0xC6,0xA9,0x7B,0x73,0x8C,0x62,0xE2,0xCB,0xAC,0x18,0x59,0x51,0x8E,0xA6,0xDC,0x61,0x71,0xD9,0x64,0x9E,0xAD,0x92,0x8F,0xC2,0x10,0x30,0xE0,0x13,0x8C,0x36,0x68,0x08,0x82,0x3D,0x1D,0x27,0xE6,0x07,0x9B,0xE6,0xE4,0xB8,0x0D,0xDD,0x91,0xCB,0x00,0x58,0x8F,0x0F,0xCC,0x60,0xA7,0x04,0x8A,0x04,0x62,0x98,0x37,0x5E,0xFC,0xAF,0x86,0xE3,0x05,0x82,0xCF,0xB0,0x71,0x03,0x5B,0xB1,0x27,0x96,0x32,0xD8,0xDA,0x6B,0xEA,0x76,0x3A,0x04,0xF8,0x62,0xD3,0xD7,0x02,0x62,0x78,0xD4,0xAA,0xEC,0x94,0x9D,0x00,0xC7,0xCF,0xB9,0xEB,0xEC,0x7B,0x19,0x84,0xB7,0x36,0x3F,0xE2,0x79,0x1A,0x51,0x18,0x8C,0xD1,0x44,0x51,0x05,0x2C,0x95,0x9D,0x78,0xEC,0x75,0x1F,0xAE,0xE1,0x3C,0x0C,0x32,0x52,0xE7,0xD9,0xDB,0xD4,0x54,0x1C,0xF0,0xA9,0x38,0x6E,0x46,0xBC,0x39,0x02,0xF0,0x4F,0x9C,0x02,0xA4,0x8E,0xCD,0x49,0xC8,0xAB,0x97,0x6C,0xA2,0xBE,0x14,0xFF,0xD1,0xF8,0xB7,0x0C,0xE4,0x05,0x49,0xD4,0xA4,0x93,0x76,0x5D,0xA9,0x82,0x7E,0x3A,0x79,0x8E,0x85,0xF7,0xAE,0x0C,0x02,0x48,0x4D,0x58,0x36,0xD0,0xE6,0x07,0x8F,0x01,0xFC,0x68,0x2B,0x5B,0x47,0x1F,0x46,0x7E,0x16,0xA7,0xF9,0x2C,0x74,0x91,0x2C,0x40,0xB9,0x10,0x3E,0x1B,0xA7,0x63,0x32,0x23,0x1F,0x30,0xB6,0x39,0x76,0xEC,0xB3,0xBF,0x7E,0x80,0xC4,0xE5,0x92,0x9C,0xD9,0x0A,0xC0,0x96,0xCD,0xDC,0x81,0x53,0x96,0x02,0x16,0xDE,0xB3,0x4D,0x6B,0x2F,0x10,0x43,0xAA,0x88,0xC4,0xEA,0x28,0x33,0xC9,0xA2,0x95,0xC3,0xA7,0xF0,0x3C,0xB8,0xD2,0x80,0x59,0x2E,0x42,0x9D,0x04,0xBF,0xBA,0xEA,0x3A,0x10,0x2B,0x1E,0xDE,0xD0,0x4A,0xA3,0x41,0x3B,0x43,0xB0,0x7D,0x78,0xCC,0x02,0x69,0xC3,0xB9,0x6A,0x00,0xAE,0x4D,0xBD,0xA5,0x9C,0x5A,0x02,0x4C,0x07,0xC8,0xE4,0x2B,0x7E,0x30,0x84,0xCF,0xD6,0xD2,0x73,0x1C,0xB8,0x44,0xBC,0x3C,0xAC,0xC9,0x99,0x44,0x85,0xCC,0x09,0x60,0xD2,0xDE,0x19,0x09,0x9E,0x24,0x5F,0x8F,0x17,0xDE,0x03,0x33,0xCD,0xAF,0x6B,0xD1,0x3C,0xEC,0x33,0x8C,0x09,0x73,0xD1,0x84,0xFE,0x98,0xD9,0x2B,0x5A,0x5F,0xCA,0x65,0xC8,0xF4,0x3C,0x40,0x7E,0x9F,0x41,0x18,0xC3,0x8E,0xC0,0x31,0x26,0x99,0x95,0x63,0x48,0x62,0x4D,0xBD,0xF9,0x99,0xB6,0xA4,0xE9,0xEB,0x01,0x00,0x77,0xE2,0xB2,0x9B,0x80,0xAF,0xA3,0x7C,0x20,0x95,0x2E,0xBC,0xE7,0x6F,0x83,0xF9,0x33,0x4B,0x3B,0xFC,0xC8,0xD9,0x62,0xB3,0x47,0xB3,0xF9,0xD2,0xEF,0x2E,0xF6,0x38,0x32,0x74,0xCB,0x2F,0x32,0xD1,0x70,0x9D,0xD2,0x6F,0x8A,0xCE,0x07,0xB1,0x01};
const unsigned char sm_lvl1[] = {0x8C,0x62,0xC7,0x9D,0x6D,0xA4,0x64,0xCD,0xC9,0xE0,0x67,0x5E,0x00,0x4F,0xCC,0x2B,0x21,0xBC,0x53,0x5E,0xA1,0xFE,0x1B,0xFA,0x0E,0x00,0xF8,0x4B,0x84,0xFF,0x75,0x98,0xBA,0xD2,0x92,0xDA,0x8B,0xCE,0x01,0xF3,0xB1,0xB9,0xF4,0xD4,0x66,0x4C,0x2F,0x4C,0x1F,0x51,0x95,0x01,0x1A,0x66,0x8A,0x32,0xE1,0x9A,0xAC,0xC8,0x0F,0xC5,0x93,0x87,0x01,0xE3,0x61,0xF3,0x85,0x32,0x57,0xB0,0x35,0x42,0xDB,0x8D,0x46,0x01,0xFD,0x23,0x8A,0xC8,0x60,0x25,0xC0,0x04,0xCA,0xEA,0x63,0x86,0x01,0x74,0x2B,0xBA,0x15,0xBE,0x2E,0x79,0x55,0x00,0x89,0x22,0xA8,0x00,0xED,0x80,0x3A,0x4A,0x7F,0xC9,0x84,0xA9,0x31,0x79,0x92,0x4A,0x01,0x49,0x0B,0x93,0xD4,0x46,0x34,0x6F,0x83,0x4C,0x5E,0x5F,0x14,0x00,0x78,0x93,0x77,0x0A,0x98,0x58,0x36,0x0F,0xB6,0xC3,0xD3,0xDF,0x00,0x3A,0x90,0x27,0xB7,0x90,0x38,0xD6,0x41,0x86,0x4F,0xFD,0xCB,0x00,0x67,0xB4,0x59,0x9A,0x7C,0x6C,0xEB,0x72,0x2F,0x90,0x6C,0xF0,0x01,0x83,0xB7,0x30,0xF4,0xF9,0x18,0x47,0xEB,0xB8,0x29,0x66,0xC6,0x01,0x5A,0x77,0xD9,0x9D,0xE4,0x56,0x7F,0x90,0x85,0x31,0x5F,0xDC,0x00,0x4E,0x27,0x43,0x58,0xDC,0x62,0xA8,0xA5,0x46,0x99,0x93,0x2A,0x00,0x00,0x95,0xF1,0x1F,0x5A,0x53,0xF3,0x53,0xB6,0x53,0x53,0x25,0x43,0xB1,0xE3,0x69,0x8F,0x84,0x90,0x5D,0x17,0x8E,0x6E,0x82,0xF2,0xD8,0x0C,0x03,0xBA,0x05,0x95,0xB0,0x9C,0xBA,0xCD,0xAF,0x0B,0xE1,0xC9,0xDA,0x00,0x9F,0x9A,0xF9,0x4E,0x99,0xCD,0xA9,0xFA,0x74,0x1E,0x64,0x40,0x34,0x00,0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
const unsigned char msg_lvl1[] = {0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
const size_t mlen_lvl1 = 33;
#elif CRYPTO_PUBLICKEYBYTES == 128
const unsigned char pk_lvl1[] = {0xCD,0x4B,0xD9,0xD6,0x2C,0xDE,0xC2,0x69,0x44,0x78,0x49,0x8E,0x3F,0x9C,0x20,0x34,0x0E,0x0C,0x36,0x43,0x71,0x8E,0x4D,0xEC,0x22,0x4A,0xCE,0x33,0x5A,0x83,0xDE,0x46,0x8A,0x27,0xC4,0x3D,0x54,0x0B,0x07,0x61,0xFE,0x13,0x71,0x7F,0xF1,0x43,0x22,0x7A,0x9F,0x48,0x49,0x3C,0x2E,0x53,0xE0,0x0F,0x18,0x8E,0xCB,0xEF,0x3E,0x39,0x12,0x00,0xF9,0x06,0x20,0xBB,0x81,0xA7,0x6D,0x8D,0x57,0xAA,0xA1,0x3B,0xDB,0xBE,0x5C,0x16,0x89,0xD7,0x66,0x83,0xDD,0x68,0xAC,0xE8,0xC9,0xC7,0xB4,0xCD,0x61,0xF3,0xAE,0x4F,0xC4,0x0A,0x87,0x94,0x9F,0xCA,0xB5,0xB7,0xBF,0xCB,0x18,0xA6,0xC7,0xCC,0xBE,0x26,0xE3,0xCA,0x67,0x88,0x02,0x89,0xA5,0x4C,0x48,0xE9,0xD8,0x42,0xDD,0x1C,0x0D,0x00};
const unsigned char sk_lvl1[] = {0xCD,0x4B,0xD9,0xD6,0x2C,0xDE,0xC2,0x69,0x44,0x78,0x49,0x8E,0x3F,0x9C,0x20,0x34,0x0E,0x0C,0x36,0x43,0x71,0x8E,0x4D,0xEC,0x22,0x4A,0xCE,0x33,0x5A,0x83,0xDE,0x46,0x8A,0x27,0xC4,0x3D,0x54,0x0B,0x07,0x61,0xFE,0x13,0x71,0x7F,0xF1,0x43,0x22,0x7A,0x9F,0x48,0x49,0x3C,0x2E,0x53,0xE0,0x0F,0x18,0x8E,0xCB,0xEF,0x3E,0x39,0x12,0x00,0xF9,0x06,0x20,0xBB,0x81,0xA7,0x6D,0x8D,0x57,0xAA,0xA1,0x3B,0xDB,0xBE,0x5C,0x16,0x89,0xD7,0x66,0x83,0xDD,0x68,0xAC,0xE8,0xC9,0xC7,0xB4,0xCD,0x61,0xF3,0xAE,0x4F,0xC4,0x0A,0x87,0x94,0x9F,0xCA,0xB5,0xB7,0xBF,0xCB,0x18,0xA6,0xC7,0xCC,0xBE,0x26,0xE3,0xCA,0x67,0x88,0x02,0x89,0xA5,0x4C,0x48,0xE9,0xD8,0x42,0xDD,0x1C,0x0D,0x00,0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x1C,0xCB,0x52,0xB5,0x61,0x25,0x56,0xF9,0x2A,0x5E,0x3D,0x5D,0xBF,0x4D,0xF3,0x3E,0xA4,0xB3,0x78,0x5A,0x6A,0xA6,0xF7,0xDA,0xF7,0x9F,0x5C,0x5F,0xC8,0x63,0x27,0x5D,0x2F,0xBC,0xFE,0x4A,0x88,0x32,0x3B,0x0F,0x24,0xF8,0x48,0x0C,0x93,0x8F,0xA7,0xDD,0x40,0x3F,0xDA,0x74,0xF7,0x86,0x4F,0x3A,0x3E,0xE6,0xF7,0x87,0x05,0x12,0x78,0xB5,0xFD,0x7E,0x2E,0xF1,0xA9,0xC4,0xD6,0x4C,0x8B,0x5F,0x8A,0x9A,0x99,0x09,0x6A,0x1E,0xF1,0x68,0x32,0x6C,0xE8,0xE1,0x12,0x21,0x44,0x05,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xC3,0x39,0xA9,0xC8,0x3E,0xB9,0x6D,0x9C,0xC2,0xD7,0x7A,0xE8,0x33,0xB3,0x01,0x26,0x7C,0x3E,0x8E,0x59,0x2C,0x19,0x65,0x92,0x35,0x84,0xF9,0x50,0x04,0x95,0x2E,0x69,0xD6,0x34,0x30,0x5E,0x4F,0x50,0xCE,0x7F,0xBD,0x89,0xE4,0x93,0x18,0x98,0x6B,0x0F,0x1F,0xA1,0xDB,0xBA,0xE1,0x7C,0x59,0x22,0x29,0x87,0x00,0x8A,0x1F,0x67,0x94,0x27,0xDA,0x64,0x5B,0xFC,0xE9,0x00,0x6C,0x68,0x15,0xE1,0x90,0x38,0xE6,0x10,0x3B,0x81,0x13,0xC2,0x58,0x7D,0xD4,0x21,0x15,0x32,0x91,0xE8,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xDB,0xE8,0xB4,0x7E,0xF0,0x4B,0x4E,0x0C,0x60,0x64,0x91,0x4B,0x85,0x8D,0xAC,0x47,0x56,0xAA,0x2C,0x89,0xD6,0x28,0xFD,0x56,0xAA,0x84,0xDD,0xBB,0xDD,0x33,0xB7,0xA8,0x48,0x5E,0xAA,0x43,0x6B,0x20,0x1A,0xFC,0x1B,0x9F,0x8A,0xEB,0xDE,0x42,0xD3,0xD8,0x96,0x20,0xB0,0x3D,0xC4,0x30,0x81,0x91,0xED,0x41,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x98,0x9F,0xA7,0xFA,0xCE,0xC9,0x8A,0xC6,0xFB,0x9E,0x4E,0x6A,0x3B,0x49,0xC4,0x4E,0xB8,0x29,0x29,0x61,0x5C,0xBA,0x18,0xD3,0xEF,0x33,0x51,0xE0,0xDA,0x72,0x0E,0xCB,0x37,0x6E,0xE1,0x4A,0x5E,0x86,0x88,0xB5,0x76,0x5E,0x24,0x1D,0xD1,0xC0,0x8A,0xC7,0x60,0x8D,0xCA,0x19,0xFB,0xB2,0x6B,0x58,0x9C,0x42,0x02,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x61,0x50,0xEB,0x78,0xD0,0xC5,0x34,0x98,0xD9,0x87,0x7C,0x87,0x6F,0x60,0x76,0xDF,0x06,0xB0,0xCB,0x7B,0x1E,0x71,0x59,0x79,0x47,0xBE,0xD1,0xF3,0x84,0x1E,0x6B,0xC1,0x47,0xEC,0x4F,0xAD,0xBA,0xA0,0xC3,0xC6,0xA5,0x8D,0x20,0xDC,0xCE,0x86,0xFC,0x21,0x80,0x91,0x25,0xC6,0xCE,0x91,0x46,0x71,0xC0,0x36,0x21,0xAA,0xD9,0xE0,0x20,0x00,0x8F,0xDD,0x84,0x7A,0x28,0x53,0x0C,0xFF,0xB9,0x6D,0xE9,0x4D,0x95,0xA2,0x28,0x4D,0x60,0x57,0xB0,0x7B,0x98,0x9C,0x0F,0x79,0xF6,0xC4,0xE4,0xDE,0x34,0xF7,0xAA,0xBD,0x44,0xF6,0xB0,0x39,0xA7,0x91,0xBB,0xD3,0x0C,0xD5,0xA7,0x73,0x16,0xB0,0xD4,0x56,0x4E,0xAD,0x4C,0x7F,0xFD,0xF7,0x86,0x0A,0x34,0x3F,0x89,0x9D,0x2E,0x1B,0x06,0x00,0x15,0xA2,0xC3,0xF1,0x82,0x15,0x53,0x70,0x3C,0xC1,0x74,0x38,0x75,0xAB,0xD4,0x98,0x4F,0xCC,0xD4,0x27,0x3C,0xD8,0x62,0xAE,0xE9,0x23,0xE7,0x7E,0x4E,0x81,0x4C,0xFF,0x5A,0xB8,0x4C,0x68,0x66,0x43,0x8A,0xF3,0xA9,0x99,0x03,0xF6,0xFA,0x21,0x7C,0x9F,0x39,0x56,0x63,0xFC,0xA3,0xB3,0x1F,0x3C,0x73,0xB4,0xCE,0x5A,0xA2,0xD2,0x18,0x00,0xE1,0xBC,0xF5,0xE6,0xDD,0xA1,0x4C,0x6D,0x6A,0x1F,0x7C,0x0D,0x1D,0xBC,0x87,0x7A,0xBD,0xBD,0x61,0xD7,0xF5,0x69,0xC7,0x77,0x8B,0xE9,0xE9,0xD7,0x5F,0x07,0x66,0x5D,0xA7,0x99,0x3E,0x01,0xE2,0xB8,0x0E,0x17,0x9F,0xA1,0x2E,0xD8,0x22,0xD4,0xBC,0x61,0xC8,0xE2,0xAB,0xF1,0x17,0x21,0x00,0xD1,0xA6,0x1E,0x3A,0x54,0x42,0x32,0x1E,0x00,0x23,0x08,0xF7,0x1A,0x86,0xE0,0x83,0x4E,0x15,0xAE,0x07,0x58,0x76,0x6E,0xAC,0x54,0x79,0x5D,0x17,0x2D,0x47,0x9C,0x6A,0x2E,0x41,0xE4,0x35,0x0C,0x59,0xDD,0xDD,0xB0,0xD0,0x7D,0xA0,0x40,0xF6,0x99,0xFA,0x22,0x2B,0xF4,0x7E,0x46,0xCF,0xE6,0x42,0x14,0x48,0xAE,0xB1,0x3E,0xE4,0xCB,0x76,0x96,0x47,0xB1,0x07,0x6E,0x7A,0x0B,0x24,0x00,0x8D,0x33,0x2C,0x3D,0x6D,0x30,0xD0,0x63,0x4C,0x19,0x48,0x2F,0xEB,0xB0,0x7D,0x8B,0x31,0xA0,0xDA,0xCA,0x31,0x0E,0x77,0x8C,0xCA,0xF4,0x01,0x89,0x24,0x8A,0xF5,0x2F,0x44,0xDC,0xC0,0x62,0x9A,0xD3,0x87,0xA0,0x56,0x2C,0x07,0xD7,0xB8,0x69,0xF6,0x2E,0x4F,0x8B,0xF3,0xAC,0x33,0xA3,0x0F,0x6E,0x9A,0x86,0xAF,0x63,0xE3,0x9D,0x1F,0x00,0xD8,0x05,0xE9,0x0F,0xB2,0x93,0x2C,0xD1,0x2B,0x90,0x05,0x28,0xD5,0x47,0xDD,0x50,0x33,0xA5,0x49,0xD4,0x26,0xF9,0x85,0xEE,0xD9,0x0D,0xFF,0x26,0x87,0xBD,0x22,0x39,0x2E,0x37,0x40,0xEB,0x86,0xD4,0x94,0x51,0xBD,0xBB,0x24,0x7B,0x4D,0x6A,0xA0,0xE8,0xEB,0xAD,0x6C,0x87,0x25,0xBB,0x10,0x58,0x94,0x53,0x60,0x13,0xA1,0x2A,0x20,0x00,0x62,0x9C,0x73,0xA7,0xAC,0x3B,0x72,0xE2,0x5F,0xA2,0xA7,0x9E,0xC8,0x15,0x96,0x86,0x92,0x40,0xD6,0x1A,0x93,0xAE,0x3B,0xE6,0x07,0x62,0xBE,0xA7,0xBD,0x5F,0x6C,0x52,0xBE,0xD7,0xC1,0x4D,0xD7,0x19,0x61,0x30,0xFC,0x71,0x52,0xFC,0x6D,0x29,0x79,0x97,0xA4,0xBB,0x6F,0x29,0xAA,0xE3,0x3E,0x63,0xB5,0xD2,0xAA,0x3A,0xF7,0x8F,0x02,0x00,0x1A,0x31,0x19,0xA5,0x55,0x27,0xC8,0x27,0x2E,0x35,0x36,0x3A,0x6C,0x1D,0xB7,0x4D,0x40,0xA3,0xFA,0x36,0x79,0xD6,0xA4,0x22,0xA8,0x9C,0x6B,0x76,0xB9,0xAA,0xC0,0x0D,0xA3,0x18,0xAB,0xD3,0x46,0xFC,0x06,0x97,0x32,0x7F,0xBC,0x8F,0x5E,0x04,0x08,0xBC,0xE1,0x5F,0x43,0x57,0x9C,0xEB,0x08,0xDE,0x03,0xF6,0x2A,0x8C,0x18,0x6D,0x24,0x00,0xE7,0xA6,0x4B,0x39,0x06,0xAD,0xE7,0x0A,0xFE,0x6B,0x12,0x40,0x2A,0x89,0xA2,0xA2,0x69,0xE6,0x60,0x9F,0x9B,0x7A,0xEC,0xAE,0x7E,0x87,0x10,0xF7,0x28,0x76,0x5E,0xB1,0xA8,0x3B,0xAA,0x5E,0x76,0xDD,0xEE,0x51,0x44,0xBC,0xF2,0xDD,0xD2,0xE3,0x45,0xB6,0x2B,0x1E,0x8F,0xD7,0xE5,0xC8,0x8B,0xFA,0x79,0xC6,0x55,0x60,0x73,0x3B,0x11,0x00,0x05,0x97,0x98,0x51,0x88,0x75,0x30,0xD1,0x5A,0xBB,0xEA,0xF0,0xE2,0x25,0x39,0x21,0x6D,0xA6,0xE3,0x36,0x44,0x23,0x8F,0xF4,0xFC,0x10,0xA8,0x34,0xA0,0xC5,0x7F,0xC9,0x4F,0xDC,0x01,0x42,0x1F,0xEE,0x2A,0x88,0x02,0x87,0x2B,0x17,0x93,0x6F,0x15,0x69,0x1B,0xC1,0x69,0x7A,0x13,0x3A,0x27,0xA4,0xC2,0x27,0x4B,0xEC,0xF7,0x7B,0x08,0x00,0x8D,0xE2,0x35,0xFC,0xD9,0xB8,0xF2,0x6D,0x53,0x41,0xC7,0x7B,0x41,0x43,0xD0,0x35,0x05,0x95,0xAC,0x09,0xA4,0x70,0x58,0x8A,0x37,0x1A,0xA4,0x57,0xBA,0x8D,0xB4,0x67,0x77,0xCB,0x8E,0x3A,0x93,0x00,0x76,0x54,0x7A,0xDD,0xA2,0xCF,0x62,0x49,0x64,0xCA,0x33,0x83,0x94,0x8A,0xD9,0x34,0x4D,0xF4,0x5E,0xD5,0xF7,0xB8,0xC2,0x85,0x0A,0x00,0x53,0x7D,0x3A,0x55,0x64,0x72,0x52,0x01,0x12,0xCC,0x09,0xE1,0x59,0x82,0x80,0x77,0x67,0x69,0x76,0xD9,0x75,0x86,0xBF,0x7E,0x99,0x45,0x0D,0x9D,0x8B,0xD9,0x37,0x18,0xDE,0x40,0xD8,0x48,0x10,0x20,0x4B,0x4B,0xA3,0x9C,0xB1,0xB4,0x33,0x55,0x12,0x16,0x07,0x27,0x18,0xB7,0x98,0xF3,0x6E,0x0A,0x38,0x06,0x07,0x49,0xA3,0x75,0x00,0x00,0x73,0xF4,0x75,0xB8,0x90,0x8A,0x7B,0x50,0x7F,0x0C,0x64,0x20,0xFD,0x45,0x78,0x05,0x09,0xE8,0x9B,0x41,0x84,0x44,0x80,0x57,0x90,0x1C,0xB4,0x06,0x42,0x9B,0xBE,0xBB,0x23,0xD5,0x42,0x96,0xCE,0xA5,0x9D,0x41,0x02,0x66,0xC5,0xE2,0x44,0x18,0xDB,0xE8,0x76,0x75,0x4A,0x96,0x09,0x40,0x3D,0x41,0x74,0x6A,0x6F,0x73,0x52,0xAC,0x0C,0x00};
const unsigned char sm_lvl1[] = {0xC8,0x56,0x5E,0x4C,0x3B,0x1C,0x20,0xF7,0x5F,0xF7,0x4C,0xFC,0xA8,0xC9,0xED,0x51,0xDB,0x97,0x01,0x5A,0xE4,0xBD,0x18,0x5F,0x4D,0xE8,0x4D,0x6B,0x1C,0x17,0x89,0xF5,0x3F,0xDA,0x0B,0xE2,0x14,0x01,0x40,0xB3,0xA7,0xE5,0x00,0x21,0x6D,0xDE,0xA2,0xB8,0x26,0xA1,0x14,0xEA,0xB8,0xC7,0x03,0x9B,0x01,0x35,0xB2,0x58,0xCE,0x60,0xCE,0xBC,0x0F,0x3D,0x40,0x40,0x0C,0x5C,0xCF,0xF6,0xCB,0xFE,0xF5,0x01,0x3D,0x56,0x62,0x1B,0xA1,0xFB,0x87,0x22,0x2D,0xD9,0xD9,0x4D,0xA3,0xE8,0xC4,0x91,0x23,0xEA,0x00,0xC9,0x21,0xC8,0x26,0xAA,0x35,0xAA,0x18,0x06,0x00,0x7B,0x49,0x27,0xB9,0x51,0xB2,0x20,0x9E,0x01,0x3A,0xF3,0xEA,0xE5,0x48,0x40,0x71,0xBE,0x33,0x77,0xE4,0xEA,0xAD,0xA2,0x08,0x53,0x1C,0x3A,0x01,0xC0,0x1F,0x64,0x3D,0x04,0x44,0xEC,0x3F,0x7A,0xF7,0xB5,0xB1,0x05,0x91,0x27,0xC6,0xE9,0x75,0x01,0xC1,0xEF,0xBD,0x52,0xA1,0x54,0xE7,0xFB,0x95,0x1F,0x93,0xCB,0xE3,0x7B,0x6C,0x1C,0x19,0xE5,0x00,0xED,0x11,0x09,0x66,0xC9,0xE0,0xAB,0x8D,0x22,0xB2,0x32,0x96,0xF5,0x11,0x28,0xB9,0x80,0x47,0x01,0x77,0x52,0xFA,0x39,0x36,0x16,0x16,0xE3,0xA2,0xD0,0x61,0x3A,0x22,0xA8,0x92,0x46,0x85,0xC6,0x01,0x77,0x87,0xAA,0x70,0x24,0x9D,0x97,0x18,0xB3,0x76,0x27,0x30,0x00,0x5B,0x36,0x4B,0x35,0x90,0x01,0x63,0xD2,0x3B,0xDB,0x01,0xDD,0x28,0x15,0xB1,0x96,0x07,0x6B,0x77,0x3D,0x1E,0x97,0xCC,0x33,0x00,0x2D,0xFD,0x24,0x8C,0xF3,0xBC,0x61,0xFC,0x28,0x8F,0x69,0xC6,0x76,0x23,0x8C,0x10,0x92,0x4D,0x01,0x00,0xE9,0xFC,0x65,0xFD,0x48,0x5E,0x13,0x6B,0x18,0x5D,0x4A,0x03,0x6C,0xDA,0xFC,0x5E,0xA0,0xBF,0xFE,0x85,0xD0,0xDF,0xB9,0x8B,0xA7,0x50,0xF3,0x12,0x83,0xCC,0x60,0x2C,0x07,0x01,0x5F,0x22,0xD1,0xEE,0xF1,0x90,0xF9,0x40,0x51,0xDC,0x35,0xCA,0x31,0x55,0x5A,0x81,0x86,0xCB,0x01,0x2B,0x84,0x71,0x95,0xEF,0x83,0x39,0x73,0x59,0xA5,0x6B,0x4E,0xF1,0xAC,0x02,0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
const unsigned char msg_lvl1[] = {0xD8,0x1C,0x4D,0x8D,0x73,0x4F,0xCB,0xFB,0xEA,0xDE,0x3D,0x3F,0x8A,0x03,0x9F,0xAA,0x2A,0x2C,0x99,0x57,0xE8,0x35,0xAD,0x55,0xB2,0x2E,0x75,0xBF,0x57,0xBB,0x55,0x6A,0xC8};
const size_t mlen_lvl1 = 33;
#else
#error "Not supported"
#endif

#if defined(ENABLE_SIGN)
static int test_sqisign_keygen() {
    int res = 0;
    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);

    res = sqisign_keypair(pk, sk);
    if (res != 0) {
        res = -1;
    }

    free(pk);
    free(sk);
    return 0;
}

static int test_sqisign_sign(int mlen) {
    int res = 0;
    unsigned char *sig = calloc(CRYPTO_BYTES + mlen, 1);
    size_t smlen = CRYPTO_BYTES + mlen;

    res = sqisign_sign(sig, &smlen, msg_lvl1, mlen, sk_lvl1);
    if (res != 0) {
        res = -1;
    }

    free(sig);
    return res;
}
#endif


static int test_sqisign_open(int mlen) {
    int res = 0;
    size_t msglen = mlen;
    unsigned char *msg = calloc(mlen, 1);

    size_t smlen = CRYPTO_BYTES + mlen;

    res = sqisign_open(msg, &msglen, sm_lvl1, smlen, pk_lvl1);
    if (res != 0) {
        res = -1;
    }


    free(msg);
    return res;
}


int main(int argc, char *argv[]) {
    int rc = 0;

    if (argc <= 1) {
        printf("Usage: <0(keygen)/1(sign)/2(open)>\n");
        exit(0);
    }

    int ar = atoi(argv[1]);

    if (ar == 0) {
#if defined(ENABLE_SIGN)
        return test_sqisign_keygen();
#else
        fprintf(stderr, "Not supported: built with signature disabled\n");
        return -1;
#endif
    } else if (ar == 1) {
#if defined(ENABLE_SIGN)
        return test_sqisign_sign(mlen_lvl1);
#else
        fprintf(stderr, "Not supported: built with signature disabled\n");
        return -1;
#endif
    } else if (ar == 2) {
        return test_sqisign_open(mlen_lvl1);
    }
    
    printf("Usage: <1(keygen)/2(sign)/3(open)>\n");
    return 0;
}
