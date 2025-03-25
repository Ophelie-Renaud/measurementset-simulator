//
// Created by hmiomand on 09/03/2020.
//

#include "fifoFunction.h"

int32_t getExponentValue(float value)
{
	uint32_t *ptr_value = (uint32_t*) &value;

	return ((*ptr_value >> FP32_MANTISSA_SIZE) & 0xff) - FP32_EXPONENT_BIAS;
}

uint64_t invert_uint64(uint64_t data_in)
{
#if 1
	uint64_t inverted_data;
	uint8_t* ptr_data_in = (uint8_t*) &data_in;
	uint8_t* ptr_inverted_data = (uint8_t*) &inverted_data;

	ptr_inverted_data[0] = ptr_data_in[7];
	ptr_inverted_data[1] = ptr_data_in[6];
	ptr_inverted_data[2] = ptr_data_in[5];
	ptr_inverted_data[3] = ptr_data_in[4];
	ptr_inverted_data[4] = ptr_data_in[3];
	ptr_inverted_data[5] = ptr_data_in[2];
	ptr_inverted_data[6] = ptr_data_in[1];
	ptr_inverted_data[7] = ptr_data_in[0];

	return inverted_data;
#else
	//CAN BE REPLACED BY return _bswap(data_in);
	return _bswap(data_in);
#endif
}

uint32_t invert_uint32(uint32_t data_in)
{
#if 1
	uint32_t inverted_data;
	uint8_t* ptr_data_in = (uint8_t*) &data_in;
	uint8_t* ptr_inverted_data = (uint8_t*) &inverted_data;

	ptr_inverted_data[0] = ptr_data_in[3];
	ptr_inverted_data[1] = ptr_data_in[2];
	ptr_inverted_data[2] = ptr_data_in[1];
	ptr_inverted_data[3] = ptr_data_in[0];

	return inverted_data;
#else
	//CAN BE REPLACED BY return _bswap(data_in);
	return _bswap(data_in);
#endif
}

uint16_t invert_uint16(uint16_t data_in)
{
	#if 1
	uint16_t inverted_data;
	uint8_t* ptr_data_in = (uint8_t*) &data_in;
	uint8_t* ptr_inverted_data = (uint8_t*) &inverted_data;

	ptr_inverted_data[0] = ptr_data_in[1];
	ptr_inverted_data[1] = ptr_data_in[0];

	return inverted_data;
	#else
	return _rotwl(data_in, 8);
	#endif
}

uint32_t fifoGet_int32_old(StorageFifo* fifo, int index)
{
	int32_t return_value = 0;
	uint8_t* data_out_arr = (uint8_t*) &return_value;

	uint32_t storageByteIndex = index * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = index * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;
	uint32_t remainingBitToRetrieve = fifo->storageSize;

	uint32_t i = 0;

	while(remainingBitToRetrieve > CHAR_BIT)
	{
		data_out_arr[workingByteIndex - (i+1)] =  (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | (fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment);
		i++;
		remainingBitToRetrieve -= CHAR_BIT;
	}

	// cas spécifique dernier/seul octet
	if (((index+1)*fifo->storageSize  & 0x07) >= (remainingBitToRetrieve & 0x07))
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1);
	}
	else
	{
		data_out_arr[workingByteIndex - (i+1)] = ((fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | ((fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1)));
	}

	return return_value;
}


void fifoSet_int32_old(StorageFifo* fifo, uint32_t data_in, int index)
{
	uint8_t* data_in_arr = (uint8_t*) &data_in;

	uint32_t storageByteIndex = index * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = index * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;

	int i = 0;

	uint32_t remainingBitToStore = fifo->storageSize;

	uint32_t widowBits = min(remainingBitToStore, rightMisalignment & 0x07);

	// If widow bits
	if (widowBits != 0)
	{
		fifo->pointerToBuffer[storageByteIndex] = ((data_in_arr[workingByteIndex - 1] >> leftMisalignment) & ((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment))
												  | (fifo->pointerToBuffer[storageByteIndex] & (uint8_t) ~((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment));
		i++;
		remainingBitToStore -= widowBits;
	}

	int j=1;

	while(remainingBitToStore > CHAR_BIT)
	{
		fifo->pointerToBuffer[storageByteIndex + i++] = (data_in_arr[workingByteIndex - j] << (rightMisalignment & 0x07))
														| (data_in_arr[workingByteIndex - j - 1] >> (CHAR_BIT - (rightMisalignment & 0x07)));
		j++;
		remainingBitToStore -= CHAR_BIT;
	}

	// cas spécifique dernier octet
	if (remainingBitToStore != 0)
	{
		if (widowBits + remainingBitToStore <= CHAR_BIT)			// If orphan bits are on the same byte
		{
			fifo->pointerToBuffer[storageByteIndex + i] = ((data_in_arr[workingByteIndex - j] << (rightMisalignment  & 0x07)) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
		else												// If orphan bits are on 2 bytes
		{
			fifo->pointerToBuffer[storageByteIndex + i] = (((data_in_arr[workingByteIndex - j] << rightMisalignment) | data_in_arr[workingByteIndex - j - 1] >> leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
	}
}


uint32_t fifoGet_int32(StorageFifo* fifo, int index)
{
	uint32_t data_out = 0;

//	int32_t storage_int32_index = index * fifo->storageSize / 32;
//	int32_t storage_int16_index = index * fifo->storageSize / 16;
//	int32_t storage_int8_index = index * fifo->storageSize / 8;

	uint32_t storage_int32_index = (index * fifo->storageSize) >> 5;
	uint32_t storage_int16_index = (index * fifo->storageSize) >> 4;
	uint32_t storage_int8_index = (index * fifo->storageSize) >> 3;

//	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) % 32);
	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) & 0x0000001F);

//	int32_t nb_bit_right = max(fifo->storageSize - nb_bit_left, 0);
	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;

	// Left side. There is always a left side

	// Retrieve data from storage (it's inverted)
	uint32_t storage_space_left_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index];

	uint32_t data_to_get_left = BYTE_SWAP32(storage_space_left_inv);

	// Putting MSBs back in place
	data_to_get_left <<= (32 - nb_bit_left);

#if 0
	if (nb_bit_left > fifo->storageSize)
	{
		// Removing bits from the next data
//		data_to_get_left &= ~(0xFFFFFFFF >> (32 - (nb_bit_left - fifo->storageSize)));
		data_to_get_left &= ~(0xFFFFFFFF >> (32 + nb_bit_right));
	}
#else
	// If nb_bit_right >= 0, then no mask is needed. Applying it anyways should do nothing
//	data_to_get_left &= ~(0xFFFFFFFFULL >> (32 + nb_bit_right));
	data_to_get_left &= ~(0xFFFFFFFFULL >> (fifo->storageSize));
#endif

	data_out = data_to_get_left;

	// Right side, if needed
	if (nb_bit_right > 0)
	{
		uint32_t storage_space_right_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1];

		uint32_t data_to_get_right = BYTE_SWAP32(storage_space_right_inv);

		// Putting LSBs back in place
		data_to_get_right >>= (nb_bit_left);

		// Removing bits from the next data
		data_to_get_right &= ~(0xFFFFFFFF >> (fifo->storageSize));

		data_out |= data_to_get_right;
	}

	return data_out;
}

void fifoSet_int32(StorageFifo* fifo, uint32_t data_in, int index)
{
	#if 0

	// Variable with _inv suffixe are bytes-inverted data, can't perform shift on them, only INVERTED masking
	// Shift operations need to be performed on non-inverted data

	// get left-side int32 index
	uint32_t storage_int32_index = (index * fifo->storageSize) >> 5;


//	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) % 32);
	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) & 0x0000001F);


//	int32_t nb_bit_right = max(fifo->storageSize - nb_bit_left, 0);
	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;


	// Removing unwanted LSBs
	data_in &= ~(0xFFFFFFFFULL >> (fifo->storageSize));

	// Getting left-side bit to store, then inverting the bytes for storage
	uint64_t data_to_store_left = ((uint64_t) data_in << 32) >> (32 - nb_bit_left);
	uint64_t data_to_store_left_inv = BYTE_SWAP64(data_to_store_left);

	// for a 64bit fetch from a 32bit index
	uint64_t storage_space_left_inv = *((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index));

	// Masking bits in memory segment
	storage_space_left_inv &= BYTE_SWAP64(~((~(0xFFFFFFFFFFFFFFFFULL >> fifo->storageSize)) >> (32 - nb_bit_left)));

	// Write left segment
//	*((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index)) = storage_space_left_inv | data_to_store_left_inv;

	atomic_fetch_and(((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index)), BYTE_SWAP64(~((~(0xFFFFFFFFFFFFFFFFULL >> fifo->storageSize)) >> (32 - nb_bit_left))));
	atomic_fetch_or(((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index)), data_to_store_left_inv);


	#else
	// Variable with _inv suffixe are bytes-inverted data, can't perform shift on them, only INVERTED masking
	// Shift operations need to be performed on non-inverted data

//	int32_t storage_int32_index = index * fifo->storageSize / 32;
//	int32_t storage_int16_index = index * fifo->storageSize / 16;
//	int32_t storage_int8_index = index * fifo->storageSize / 8;

	uint32_t storage_int32_index = (index * fifo->storageSize) >> 5;
	uint32_t storage_int16_index = (index * fifo->storageSize) >> 4;
	uint32_t storage_int8_index = (index * fifo->storageSize) >> 3;

//	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) % 32);
	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) & 0x0000001F);


//	int32_t nb_bit_right = max(fifo->storageSize - nb_bit_left, 0);
	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;

	// Left side. There is always a left side

	// Removing unwanted LSBs
	data_in &= ~(0xFFFFFFFFULL >> (fifo->storageSize));

	// Getting left-side bit to store, then inverting the bytes for storage
	uint32_t data_to_store_left = data_in >> (32 - nb_bit_left);
	uint32_t data_to_store_left_inv = BYTE_SWAP32(data_to_store_left);

	#if 1
	uint32_t storage_space_left_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index];

	uint32_t storage_space_left_inv_prep = storage_space_left_inv & BYTE_SWAP32(~(0xFFFFFFFF >> (32 - nb_bit_left)) | (0xFFFFFFFFULL >> (32 - (nb_bit_left - fifo->storageSize))));

	// Write left segment
	((uint32_t*) fifo->pointerToBuffer)[storage_int32_index] = storage_space_left_inv_prep | data_to_store_left_inv;

	#else
	uint32_t storage_space_left_inv;
	uint32_t storage_space_left_inv_prep;

	do
	{

		// for a 32bit fetch
		storage_space_left_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index];

		storage_space_left_inv_prep = storage_space_left_inv & BYTE_SWAP32(~(0xFFFFFFFF >> (32 - nb_bit_left)) | (0xFFFFFFFFULL >> (32 - (nb_bit_left - fifo->storageSize))));

		// Write left segment
//		((uint32_t*) fifo->pointerToBuffer)[storage_int32_index] = storage_space_left_inv_prep | data_to_store_left_inv;

	}while(atomic_compare_exchange_strong(&((uint32_t*) fifo->pointerToBuffer)[storage_int32_index], &storage_space_left_inv, storage_space_left_inv_prep | data_to_store_left_inv) == 0);
	#endif

	// Right side, if needed
	if (nb_bit_right > 0)
	{

		// scrapping already stored MSBs
		uint32_t data_to_store_right = data_in << nb_bit_left;

		// scrapping unwanted LSBs
//			data_to_store_right &= ~((0x01 << (32 - nb_bit_right)) - 1);
		data_to_store_right &= ~(0xFFFFFFFFULL >> nb_bit_right);

		// inverting data for coming insertion
		uint32_t data_to_store_right_inv = BYTE_SWAP32(data_to_store_right);

		#if 1
		// Fetch the right side 32 bits segment
		uint32_t storage_space_right_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index + 1];

		// Cleaning destination
		storage_space_right_inv &= BYTE_SWAP32(0xFFFFFFFFULL >> nb_bit_right);

		// Write data to memory
		((uint32_t*) fifo->pointerToBuffer)[storage_int32_index + 1] =
				storage_space_right_inv | data_to_store_right_inv;

//			__sync_fetch_and_and(&((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1], BYTE_SWAP32(0xFFFFFFFFULL >> nb_bit_right));
//			__sync_fetch_and_or(&((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1], data_to_store_right_inv);

//			atomic_fetch_and(&((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1], BYTE_SWAP32(0xFFFFFFFFULL >> nb_bit_right));
//			atomic_fetch_or(&((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1], data_to_store_right_inv);
		#else
		uint32_t storage_space_right_inv;
		uint32_t storage_space_right_inv_prep;

		do
		{

			// for a 32bit fetch
			storage_space_right_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1];

			storage_space_right_inv_prep = storage_space_right_inv & BYTE_SWAP32(0xFFFFFFFFULL >> nb_bit_right);

			// Write left segment
//		((uint32_t*) fifo->pointerToBuffer)[storage_int32_index] = storage_space_left_inv_prep | data_to_store_left_inv;

		}while(atomic_compare_exchange_strong(&((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1], &storage_space_right_inv, storage_space_right_inv_prep | data_to_store_right_inv) == 0);
		#endif


		// Releasing lock ?
//		}
	}
	#endif
}

uint32_t fifoGet_int32_8bytes(StorageFifo* fifo, int index)
{
	// Seems to be only efficient for 31bit wide data, and for some reason also 32bits data (probably because it doesn't have conditional jump)

	// Variable with _inv suffixe are bytes-inverted data, can't perform shift on them, only INVERTED masking
	// Shift operations need to be performed on non-inverted data
	// Read data from AxB with a 8bytes chunk with 4bytes alignment

	uint32_t data_out = 0;

	uint32_t storage_int32_index = (index * fifo->storageSize) >> 5;
	uint32_t storage_int16_index = (index * fifo->storageSize) >> 4;
	uint32_t storage_int8_index = (index * fifo->storageSize) >> 3;

	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) & 0x0000001F);

	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;

	// Retrieve data from storage (it's inverted)
	uint64_t storage_space_left_inv = *((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index));

	uint64_t data_to_get_left = BYTE_SWAP64(storage_space_left_inv);

	// Putting MSBs back in place
	data_to_get_left = (data_to_get_left << (32 - nb_bit_left)) >> 32;

	return data_to_get_left;
}


void fifoSet_int32_8bytes(StorageFifo* fifo, uint32_t data_in, int index)
{
	// Variable with _inv suffixe are bytes-inverted data, can't perform shift on them, only INVERTED masking
	// Shift operations need to be performed on non-inverted data
	// Write data in AxB with a 8bytes chunk with 4bytes alignment

	// get left-side int32 index
	uint32_t storage_int32_index = (index * fifo->storageSize) >> 5;

	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) & 0x0000001F);

	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;

	// for a 64bit fetch from a 32bit index
	uint64_t storage_space_left_inv = *((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index));

	// Removing unwanted LSBs
	data_in &= ~(0xFFFFFFFFULL >> (fifo->storageSize));

	// Getting left-side bit to store, then inverting the bytes for storage
	uint64_t data_to_store_left = ((uint64_t) data_in << 32) >> (32 - nb_bit_left);
	uint64_t data_to_store_left_inv = BYTE_SWAP64(data_to_store_left);

	// Masking bits in memory segment
	storage_space_left_inv &= BYTE_SWAP64(~((~(0xFFFFFFFFFFFFFFFFULL >> fifo->storageSize)) >> (32 - nb_bit_left)));

	// Write left segment
	*((uint64_t*) (((uint32_t*) (fifo->pointerToBuffer)) + storage_int32_index)) = storage_space_left_inv | data_to_store_left_inv;
}

void fifoGet(StorageFifo* fifo, void* data_out, int index)
{
//	char* data_out_arr = data_out;
	uint8_t* data_out_arr = (uint8_t*) data_out;

	// cleaning output just in case
	for (int i = 0; i < fifo->workingSize/CHAR_BIT; i++)
		data_out_arr[i] = 0x00;

	uint32_t storageByteIndex = index * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = index * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;
	uint32_t remainingBitToRetrieve = fifo->storageSize;

	uint32_t i = 0;

	while(remainingBitToRetrieve > CHAR_BIT)
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												 | (fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment);
		i++;
		remainingBitToRetrieve -= CHAR_BIT;
	}

	// cas spécifique dernier/seul octet
	if (((index+1)*fifo->storageSize  & 0x07) >= (remainingBitToRetrieve & 0x07))
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1);
	}
	else
	{
		data_out_arr[workingByteIndex - (i+1)] = ((fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | ((fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1)));
	}
}

float fifoGet_float(StorageFifo* fifo, int index)
{
	float return_value;

	// Adjust index in case of multi-instance AxB
	index += fifo->indexOffset;

	switch (fifo->fifoType)
	{
		case FxP :
			return_value = fifoGet_float_FxP(fifo, index);
			break;
		case FxP_control :
			return_value = fifoGet_float_FxP(fifo, index);
			break;
		case cFP :
			return_value = fifoGet_float_cFP(fifo, index);
			break;
		case SCALING :
			return_value = fifoGet_float_scaling(fifo, index);
			break;
		case SCALING_control :
			return_value = fifoGet_float_scaling(fifo, index);
			break;
		case NO_TYPE :
		{
			uint32_t temp = FIFO_GET_FUNC(fifo, index);
//			uint32_t temp = fifoGet_int32_16bits(fifo, index);
//			uint32_t temp = fifo->ptr_fifoGet_int32(fifo, index);
			return_value = *((float*) &temp);
			break;
		}
		default :
		{
			uint32_t temp = FIFO_GET_FUNC(fifo, index);
//			uint32_t temp = fifo->ptr_fifoGet_int32(fifo, index);
			return_value = *((float*) &temp);
		}
	}


	return return_value;
}

float fifoSet_float(StorageFifo* fifo, float data_in, int index)
{

	// Adjust index in case of multi-instance AxB
	index += fifo->indexOffset;

	switch (fifo->fifoType)
	{
		case FxP :
			fifoSet_float_FxP(fifo, data_in, index);
			break;
		case FxP_control :
			fifoSet_float_FxP_control(fifo, data_in, index);
			break;
		case cFP :
			fifoSet_float_cFP(fifo, data_in, index);
			break;
		case SCALING :
			fifoSet_float_scaling(fifo, data_in, index);
//			fifoSet_float_scaling_control(fifo, data_in, index);
			break;
		case SCALING_control :
//			fifoSet_float_scaling(fifo, data_in, index);
			fifoSet_float_scaling_control(fifo, data_in, index);
			break;
		case NO_TYPE :
		{
			uint32_t *temp = (uint32_t*) &data_in;
			FIFO_SET_FUNC(fifo, *temp, index);
			break;
		}
		default :
		{
			uint32_t *temp = (uint32_t*) &data_in;
			FIFO_SET_FUNC(fifo, *temp, index);
		}
	}

	return data_in;
}

void fifoSet(StorageFifo* fifo, void* data_in, int index)
{
	uint8_t* data_in_arr = (uint8_t*) data_in;

	uint32_t storageByteIndex = index * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = index * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;

	int i = 0;

	uint32_t remainingBitToStore = fifo->storageSize;

	uint32_t widowBits = min(remainingBitToStore, rightMisalignment & 0x07);

	// If widow bits
	if (widowBits != 0)
	{
		fifo->pointerToBuffer[storageByteIndex] = ((data_in_arr[workingByteIndex - 1] >> leftMisalignment) & ((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment))
												  | (fifo->pointerToBuffer[storageByteIndex] & (uint8_t) ~((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment));
		i++;
		remainingBitToStore -= widowBits;
	}

	int j=1;

	while(remainingBitToStore > CHAR_BIT)
	{
		fifo->pointerToBuffer[storageByteIndex + i++] = (data_in_arr[workingByteIndex - j] << (rightMisalignment & 0x07))
														| (data_in_arr[workingByteIndex - j - 1] >> (CHAR_BIT - (rightMisalignment & 0x07)));
		j++;
		remainingBitToStore -= CHAR_BIT;
	}

	// cas spécifique dernier octet
	if (remainingBitToStore != 0)
	{
		if (widowBits + remainingBitToStore <= CHAR_BIT)			// If orphan bits are on the same byte
		{
			fifo->pointerToBuffer[storageByteIndex + i] = ((data_in_arr[workingByteIndex - j] << (rightMisalignment  & 0x07)) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
		else												// If orphan bits are on 2 bytes
		{
			fifo->pointerToBuffer[storageByteIndex + i] = (((data_in_arr[workingByteIndex - j] << rightMisalignment) | data_in_arr[workingByteIndex - j - 1] >> leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
	}
}

void fifoGet_override(StorageFifo* fifo, void* data_out, int index, uint32_t storageSize, uint32_t workingSize)
{
//	char* data_out_arr = data_out;
	uint8_t* data_out_arr = (uint8_t*) data_out;

	// cleaning output just in case
	for (int i = 0; i < workingSize/CHAR_BIT; i++)
		data_out_arr[i] = 0x00;

	uint32_t storageByteIndex = index * storageSize >> 3;
	uint32_t workingByteIndex = workingSize >> 3;
	uint8_t leftMisalignment = index * storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;
	uint32_t remainingBitToRetrieve = storageSize;

	uint32_t i = 0;

	while(remainingBitToRetrieve > CHAR_BIT)
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												 | (fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment);
		i++;
		remainingBitToRetrieve -= CHAR_BIT;
	}

	// cas spécifique dernier/seul octet
	if (((index+1)*storageSize  & 0x07) >= (remainingBitToRetrieve & 0x07))
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1);
	}
	else
	{
		data_out_arr[workingByteIndex - (i+1)] = ((fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | ((fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1)));
	}
}

void fifoSet_override(StorageFifo* fifo, void* data_in, int index, uint32_t storageSize, uint32_t workingSize)
{
	uint8_t* data_in_arr = (uint8_t*) data_in;

	uint32_t storageByteIndex = index * storageSize >> 3;
	uint32_t workingByteIndex = workingSize >> 3;
	uint8_t leftMisalignment = index * storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;

	int i = 0;

	uint32_t remainingBitToStore = storageSize;

	uint32_t widowBits = min(remainingBitToStore, rightMisalignment & 0x07);

	// If widow bits
	if (widowBits != 0)
	{
		fifo->pointerToBuffer[storageByteIndex] = ((data_in_arr[workingByteIndex - 1] >> leftMisalignment) & ((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment))
												  | (fifo->pointerToBuffer[storageByteIndex] & (uint8_t) ~((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment));
		i++;
		remainingBitToStore -= widowBits;
	}

	int j=1;

	while(remainingBitToStore > CHAR_BIT)
	{
		fifo->pointerToBuffer[storageByteIndex + i++] = (data_in_arr[workingByteIndex - j] << (rightMisalignment & 0x07))
														| (data_in_arr[workingByteIndex - j - 1] >> (CHAR_BIT - (rightMisalignment & 0x07)));
		j++;
		remainingBitToStore -= CHAR_BIT;
	}

	// cas spécifique dernier octet
	if (remainingBitToStore != 0)
	{
		if (widowBits + remainingBitToStore <= CHAR_BIT)			// If orphan bits are on the same byte
		{
			fifo->pointerToBuffer[storageByteIndex + i] = ((data_in_arr[workingByteIndex - j] << (rightMisalignment  & 0x07)) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
		else												// If orphan bits are on 2 bytes
		{
			fifo->pointerToBuffer[storageByteIndex + i] = (((data_in_arr[workingByteIndex - j] << rightMisalignment) | data_in_arr[workingByteIndex - j - 1] >> leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
	}
}

void fifoSimulate_float(int NB_ELEMENT, float* buffer, int storageSize)
{
	float error = 0.0f;

	float temp;
	uint32_t* ptr_temp = &temp;
	uint32_t* ptr_buffer_i;

	StorageFifo fifo = fifoCreate(NB_ELEMENT, storageSize, 32);

	for (int i = 0; i < NB_ELEMENT; i++)
	{
		fifoSet_float(&fifo, buffer[i], i);
	}

	for (int i = 0; i < NB_ELEMENT; i++)
	{
		temp = fifoGet_float(&fifo, i);

		error += fabsf(temp - buffer[i]);

//		ptr_buffer_i = &buffer[i];
//		if (*ptr_temp != ((*ptr_buffer_i) & ~(0xFFFFFFFFULL >> (fifo.storageSize))))
//			__asm__("nop");

		buffer[i] = temp;
	}

	if (error != 0)
	{
		printf("Error = %f\n", error/NB_ELEMENT);
		__asm__("nop");
	}


	fifoDestroy(&fifo);
}



StorageFifo fifoCreate(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	StorageFifo fifo;

	// malloc to the upper 8 bytes
	uint64_t byteToAlloc = (nb_element * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

//	fifo.pointerToBuffer = (uint8_t*) malloc(((nb_element * storageSize / CHAR_BIT) + (64 - 1)) & ~(64 - 1));
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;

	fifo.fifoType = NO_TYPE;

//	if (fifo.storageSize > 30)
//		fifo.ptr_fifoGet_int32 = &fifoGet_int32_8bytes;
//	else
//		fifo.ptr_fifoGet_int32 = &fifoGet_int32;

	fifo.indexOffset = 0;
	fifo.nbInstance = 1;

	return fifo;
}

StorageFifo fifoCreate_multi(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, int nb_instance)
{
//	StorageFifo AxB = fifoCreate(nb_element, storageSize, workingSize);
	StorageFifo AxB = fifoCreate_FxP(nb_element, storageSize, workingSize, 120);

	AxB.nbInstance = nb_instance;

	return AxB;
}

StorageFifo fifoCreate_slave(StorageFifo* master_AxB, int actor_instance)
{
	// Copy info from master
	StorageFifo slave_AxB = *master_AxB;

	slave_AxB.indexOffset = (slave_AxB.nb_element / slave_AxB.nbInstance) * actor_instance;


	return slave_AxB;
}

void fifoDestroy(StorageFifo* fifo)
{
	if (fifo->pointerToBuffer != NULL)
	{
		free(fifo->pointerToBuffer);
		fifo->pointerToBuffer = NULL;
	}
}

void fifoReset(StorageFifo* fifo)
{
	switch (fifo->fifoType)
	{
		case FxP :
			__asm__("nop");
			break;
		case FxP_control :
			fifoReset_FxP(fifo);
			break;
		case cFP :
			__asm__("nop");
			break;
		case SCALING :
			__asm__("nop");
			break;
		case SCALING_control :
			fifoReset_scaling(fifo);
			break;
		default :
			__asm__("nop");
	}
}



// FxP

StorageFifo fifoCreate_FxP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float maxValue)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	int storageIntegerSize = (int) (floorf(log2f(maxValue)) +2);

	int storageDecimalSize = storageSize - storageIntegerSize;

	StorageFifo fifo;

	uint64_t byteToAlloc = (nb_element * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

	//fifo.pointerToBuffer = (uint8_t*) malloc(nb_element * storageSize / CHAR_BIT);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);

	fifo.storageDecimalSize = storageDecimalSize;
	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;

	fifo.originalDecimalSize = storageDecimalSize;

	fifo.fifoType = FxP;
//	fifo.fifoType = FxP_control;

	fifo.nbInstance = 1;
	fifo.indexOffset = 0;

	fifo.currentStorageIntegerSize = INT32_MIN;

	return fifo;
}

float fifoGet_float_FxP(StorageFifo* fifo, int index)
{
	float return_value = 0.0f;

	int32_t temp = FIFO_GET_FUNC(fifo, index);

	temp >>= (fifo->workingSize - fifo->storageSize);

	return_value = (float) temp;

	// It's faster to directly modify the exponent field of the float than computing the appropriate power of 2
	if (return_value != 0)
	{
		uint32_t* temp = (uint32_t*) &return_value;

		uint8_t exponent = (*temp >> FP32_MANTISSA_SIZE) & 0xff;
		exponent -= fifo->storageDecimalSize;

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << FP32_MANTISSA_SIZE);
	}

	return return_value;
}


void fifoSet_float_FxP(StorageFifo* fifo, float data_in, int index)
{
	// It's faster to directly modify the exponent field of the float than computing the appropriate power of 2
	if (data_in != 0)
	{
		uint32_t* temp = (uint32_t*) &data_in;

		uint8_t exponent = (*temp >> FP32_MANTISSA_SIZE) & 0xff;

		exponent += fifo->storageDecimalSize;

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << FP32_MANTISSA_SIZE);
	}

	// Casting data_in to int32 to complete FxP conversion
	int32_t temp = (int32_t) data_in;

	// Setting temp as right-padded to work with fifoSet_int32
	temp <<= (fifo->workingSize - fifo->storageSize);

	FIFO_SET_FUNC(fifo, temp, index);
}

void fifoSimulate_float_FxP(int NB_ELEMENT, float* buffer, int storageSize)
{
	float max_value = 0.0f;
//	float min_value = INFINITY;

	// Get highest value
	for (int i = 0; i < NB_ELEMENT; i++)
	{
		if (fabsf(buffer[i]) > max_value)
			max_value = fabsf(buffer[i]);
	}

	StorageFifo fifo = fifoCreate_FxP(NB_ELEMENT, storageSize, 32, max_value);

	for (int i = 0; i < NB_ELEMENT; i++)
		fifoSet_float(&fifo, buffer[i], i);

	for (int i = 0; i < NB_ELEMENT; i++)
		buffer[i] = fifoGet_float(&fifo, i);

	fifoDestroy(&fifo);
}


void fifoSet_float_FxP_control(StorageFifo* fifo, float data_in, int index)
{

	int32_t storageIntegerSize = ((int) (floorf(log2f(data_in)) +2 ));

	if (storageIntegerSize > fifo->currentStorageIntegerSize)
		fifo->currentStorageIntegerSize = storageIntegerSize;

	// Check if data_in can be stored with fifo->storageDecimalSize
	if (storageIntegerSize > (fifo->storageSize - fifo->storageDecimalSize))
	{
		// if data_in cannot be stored as is

		// determine new fifo->storageDecimalSize but do not change it yet
		int32_t newStorageDecimalSize = fifo->storageSize - storageIntegerSize;

		for (int i = 0; i < fifo->nb_element; i++)
		{
			float temp = fifoGet_float_FxP(fifo, i);
			fifoSet_float_FxP_override(fifo, temp, i, newStorageDecimalSize);

			float check = fifoGet_float_FxP_override(fifo, i, newStorageDecimalSize);
			__asm__("nop");
		}
		fifo->storageDecimalSize = newStorageDecimalSize;
	}

	fifoSet_float_FxP(fifo, data_in, index);
}

float fifoGet_float_FxP_override(StorageFifo* fifo, int index, uint32_t storageDecimalSize)
{
	float return_value = 0.0f;

	// fifoGet_int32 return data that is right-padded with 0s
	return_value = (float) (FIFO_GET_FUNC(fifo, index) >> (fifo->workingSize - fifo->storageSize));

	if (return_value != 0)
	{
		uint32_t* temp = (uint32_t*) &return_value;

		uint8_t exponent = (*temp >> FP32_MANTISSA_SIZE) & 0xff;
		exponent -= storageDecimalSize;

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << FP32_MANTISSA_SIZE);
	}

	__asm__("nop");

	return return_value;
}


void fifoSet_float_FxP_override(StorageFifo* fifo, float data_in, int index, uint32_t storageDecimalSize)
{
//	int32_t temp_buffer;

//	uint8_t* data_in_arr = (uint8_t*) &temp_buffer;

//	printf("%x\n", *((uint32_t*) &data_in));

	if (data_in != 0)
	{
		// Creating pointer to perform binary operation on float
		uint32_t* temp = (uint32_t*) &data_in;

		// Get exponent value
		uint8_t exponent = (*temp >> FP32_MANTISSA_SIZE) & 0xff;

		exponent += storageDecimalSize;

		// Applying exponent value
		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << FP32_MANTISSA_SIZE);
	}

	// Casting data_in to int32 to complete FxP conversion
	int32_t temp = (int32_t) data_in;

	// Setting temp as right-padded to work with fifoSet_int32
	temp <<= (fifo->workingSize - fifo->storageSize);

	FIFO_SET_FUNC(fifo, temp, index);
}

void fifoReset_FxP(StorageFifo* fifo)
{
//	fifo->storageDecimalSize = fifo->originalDecimalSize;

	if (fifo->currentStorageIntegerSize < (fifo->storageSize - fifo->storageDecimalSize))
	{
		fifo->storageDecimalSize = fifo->storageSize - fifo->currentStorageIntegerSize;
	}
}





StorageFifo fifoCreate_scaling(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float minValue, float maxValue)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	float amplitude = maxValue - minValue;

	StorageFifo fifo;

	uint64_t byteToAlloc = (nb_element * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

	//fifo.pointerToBuffer = (uint8_t*) malloc(nb_element * storageSize / CHAR_BIT);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;

	fifo.fifoType = SCALING;
//	fifo.fifoType = SCALING_control;

	fifo.offset = -minValue;
//	fifo.stepSize = amplitude * 1.01 / (powf(2, storageSize));
	fifo.stepSize = amplitude * SCALING_AMPLITUDE_CORRECTION / (powf(2, storageSize));

	fifo.offset = (int32_t) (-minValue / fifo.stepSize);


	fifo.originalStepSize = fifo.stepSize;
	fifo.originalOffset = fifo.offset;

	fifo.minValue = minValue;
	fifo.maxValue = maxValue;

	fifo.currentMaxValue = -INFINITY;
	fifo.currentMinValue = INFINITY;

	fifo.indexOffset = 0;
	fifo.nbInstance = 1;

	return fifo;
}

float fifoGet_float_scaling(StorageFifo* fifo, int index)
{
	float return_value = 0.0f;

	// fifoGet_int32 return data that is right-padded with 0s
	return_value = (float) ((uint32_t) FIFO_GET_FUNC(fifo, index) >> (fifo->workingSize - fifo->storageSize));

//	return_value = return_value*fifo->stepSize - fifo->offset;

	return_value = (return_value - fifo->offset) * fifo->stepSize;

	__asm__("nop");

	return return_value;
}


void fifoSet_float_scaling(StorageFifo* fifo, float data_in, int index)
{
	// Applying scaling
//	data_in = (data_in + fifo->offset) / fifo->stepSize;

	data_in = data_in/fifo->stepSize + fifo->offset;

	// Casting data_in to int32 to complete FxP conversion
	uint32_t temp = (uint32_t) data_in;

	// Dirty patch, need to actually fix the stepsize
	if (temp >= ((uint64_t) 0x01 << fifo->storageSize))
	{
		__asm__("nop");
//		temp = (0x01 << fifo->storageSize) -1;
		printf("Overflow detected, please raise max value\n");
	}

	// Setting temp as right-padded to work with fifoSet_int32
	temp <<= (fifo->workingSize - fifo->storageSize);

	FIFO_SET_FUNC(fifo, temp, index);
}

void fifoSimulate_float_scaling(int NB_ELEMENT, float* buffer, int storageSize)
{
	float max_value = -INFINITY;
	float min_value = INFINITY;

	// Get highest value
	for (int i = 0; i < NB_ELEMENT; i++)
	{
		if (buffer[i] > max_value)
			max_value = buffer[i];

		if (buffer[i] < min_value)
			min_value = buffer[i];
	}

	StorageFifo fifo = fifoCreate_scaling(NB_ELEMENT, storageSize, 32, min_value, max_value);

	for (int i = 0; i < NB_ELEMENT; i++)
	{
		fifoSet_float(&fifo, buffer[i], i);
	}

	for (int i = 0; i < NB_ELEMENT; i++)
	{
		buffer[i] = fifoGet_float(&fifo, i);
	}

	fifoDestroy(&fifo);
}

void fifoSet_float_scaling_control(StorageFifo* fifo, float data_in, int index)
{
//	float minValue = -fifo->offset;
//	float maxValue = fifo->stepSize*powf(2, fifo->storageSize) + minValue;

	if (data_in < fifo->currentMinValue)
		fifo->currentMinValue = data_in;
	if (data_in > fifo->currentMaxValue)
		fifo->currentMaxValue = data_in;


	// Check if current value is lower than the lowest one, if so:
	if (data_in < fifo->minValue)
	{
		__asm__("nop");

		fifo->minValue = data_in;

		// Prepare newOffset value and newStepSize value
		float newOffset = -fifo->minValue;
		float newStepSize = (fifo->maxValue - fifo->minValue) * SCALING_AMPLITUDE_CORRECTION / powf(2, fifo->storageSize);

		// Read fifo value with old parameters, store them back with new parameter
		for (int i = 0; i < fifo->nb_element; i++)
		{
			float temp = fifoGet_float_scaling(fifo, i);
			fifoSet_float_scaling_override(fifo, temp, i, newStepSize, newOffset);

//			fifoSet_float_scaling_override(fifo, fifoGet_float_scaling(fifo, i), i, newStepSize, newOffset);

			float check = fifoGet_float_scaling_override(fifo, i, newStepSize, newOffset);
			__asm__("nop");
		}

		// Set new parameters as default
		fifo->stepSize = newStepSize;
		fifo->offset = newOffset;
	}

	// Check if current value is higher than the highest one, if so:
	if (data_in >= fifo->maxValue)
	{
		__asm__("nop");

		fifo->maxValue = data_in;

		// Prepare newStepSize value. 1.01 is a small offset to prevent maxValue to be encoded the same as minValue.
		// Appropriate offset needs to be determined
		float newStepSize = (fifo->maxValue - fifo->minValue) * SCALING_AMPLITUDE_CORRECTION / (powf(2, fifo->storageSize));

		// Read fifo value with old parameters, store them back with new parameter
		for (int i = 0; i < fifo->nb_element; i++)
		{
			float temp = fifoGet_float_scaling(fifo, i);
			fifoSet_float_scaling_override(fifo, temp, i, newStepSize, fifo->offset);

//			fifoSet_float_scaling_override(fifo, fifoGet_float_scaling(fifo, i), i, newStepSize, fifo->offset);

			float check = fifoGet_float_scaling_override(fifo, i, newStepSize, fifo->offset);
			__asm__("nop");
		}

		// Set new parameter as default
		fifo->stepSize = newStepSize;
	}

	fifoSet_float_scaling(fifo, data_in, index);
}

float fifoGet_float_scaling_override(StorageFifo* fifo, int index, float stepSize, float offset)
{
	float return_value = 0.0f;

	// fifoGet_int32 return data that is right-padded with 0s
	return_value = (float) ((uint32_t) FIFO_GET_FUNC(fifo, index) >> (fifo->workingSize - fifo->storageSize));

	return_value = return_value*stepSize - offset;

	__asm__("nop");

	return return_value;
}

void fifoSet_float_scaling_override(StorageFifo* fifo, float data_in, int index, float stepSize, float offset)
{
	// Applying scaling
	data_in = (data_in + offset) / stepSize;

	// Casting data_in to int32 to complete FxP conversion
	uint32_t temp = (uint32_t) data_in;

	// Dirty patch, need to actually fix the stepsize
	if (temp >= ((uint64_t) 0x01 << (fifo->storageSize)))
	{
		__asm__("nop");
//		temp = (0x01 << fifo->storageSize) -1;
		printf("Overflow detected, please raise max value\n");
	}

	// Setting temp as right-padded to work with fifoSet_int32
	temp <<= (fifo->workingSize - fifo->storageSize);

	FIFO_SET_FUNC(fifo, temp, index);
}

void fifoReset_scaling(StorageFifo* fifo)
{
//	fifo->offset = fifo->originalOffset;
//	fifo->stepSize = fifo->originalStepSize;

	if (fifo->fifoType == SCALING_control)
	{
		float currentAmplitude = (fifo->currentMaxValue - fifo->currentMinValue);
		float margin = currentAmplitude * 0.05f;

		if ((fifo->maxValue - fifo->minValue) > (currentAmplitude + margin))
		{
//			fifo->minValue = fifo->currentMinValue - margin / 2;
//			fifo->maxValue = fifo->currentMaxValue + margin / 2;
			fifo->minValue = fifo->currentMinValue;
			fifo->maxValue = fifo->currentMaxValue;

			fifo->offset = -fifo->minValue;
			fifo->stepSize = (fifo->maxValue - fifo->minValue) * SCALING_AMPLITUDE_CORRECTION / (powf(2, fifo->storageSize));
		}

		fifo->currentMaxValue = -INFINITY;
		fifo->currentMinValue = INFINITY;
	}
}


// cFP

StorageFifo fifoCreate_cFP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float maxAbsExponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	int exponentSize = (int) (floorf(log2f(maxAbsExponent)+1)) +1;
	uint32_t exponentBias = (uint32_t) powf(2,exponentSize-1)-1;
	int mantissaSize = storageSize - exponentSize -1;

	if (mantissaSize > FP32_MANTISSA_SIZE)
		mantissaSize = FP32_MANTISSA_SIZE;

	StorageFifo fifo;

	fifo.pointerToBuffer = (uint8_t*) malloc(nb_element * storageSize / CHAR_BIT);

	fifo.exponentBias = exponentBias;
	fifo.exponentSize = exponentSize;
	fifo.mantissaSize = mantissaSize;
	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.signBit = SIGN_UNDEF;
	fifo.nb_element = nb_element;

	fifo.fifoType = cFP;

	return fifo;
}

StorageFifo fifoCreate_cFP_better(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float minExponent, float maxExponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}



	// +1 is extra room to encode 0
	int exponentSize = (int) (floorf(log2f(maxExponent - minExponent +1)+1));
	// +1 is extra room to encode 0
	int32_t exponentBias = -minExponent +1;

	int mantissaSize = storageSize - exponentSize -1;

	if (mantissaSize > FP32_MANTISSA_SIZE)
		mantissaSize = FP32_MANTISSA_SIZE;

	StorageFifo fifo;

	uint64_t byteToAlloc = (nb_element * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

//	fifo.pointerToBuffer = (uint8_t*) malloc(nb_element * storageSize / CHAR_BIT);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);

	fifo.exponentBias = exponentBias;
	fifo.exponentSize = exponentSize;
	fifo.mantissaSize = mantissaSize;
	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.signBit = SIGN_UNDEF;

	fifo.fifoType = cFP;

	fifo.indexOffset = 0;
	fifo.nbInstance = 1;

	return fifo;
}

StorageFifo fifoCreate_cFP_better_sign(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float minExponent, float maxExponent, int sign)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}



	// +1 is extra room to encode 0, might need to add +1 to encode inf and NaN
	int exponentSize = (int) (floorf(log2f(maxExponent - minExponent +1)+1));
	// +1 is extra room to encode 0
	int32_t exponentBias = -minExponent +1;

	int mantissaSize = storageSize - exponentSize;

	// If buffer sign is not uniform, add a sign bit
	if (sign == SIGN_UNDEF)
		mantissaSize--;

	if (mantissaSize > FP32_MANTISSA_SIZE)
		mantissaSize = FP32_MANTISSA_SIZE;

	StorageFifo fifo;

	fifo.pointerToBuffer = (uint8_t*) malloc(nb_element * storageSize / CHAR_BIT);

	fifo.exponentBias = exponentBias;
	fifo.exponentSize = exponentSize;
	fifo.mantissaSize = mantissaSize;
	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.signBit = sign;

	fifo.fifoType = cFP;

	return fifo;
}

float fifoGet_float_cFP(StorageFifo* fifo, int index)
{
	float return_value_cFP = 0.0f;

	// fifoGet_int32 return data that is right-padded with 0s
	int32_t* return_value = (int32_t*) &return_value_cFP;

	*return_value = FIFO_GET_FUNC(fifo, index);

	if (return_value_cFP != 0)
	{
		// Get pointer
		uint32_t* temp = (uint32_t*) &return_value_cFP;

		// Get sign
		uint32_t sign = *temp & 0x80000000;

		// Get exponent, on 8 bits or less
		int32_t exponent = (*temp >> (fifo->mantissaSize + (fifo->workingSize - fifo->storageSize)));

		// Isolate exponent bits
		exponent &= ((1<<(fifo->exponentSize))-1);

		// Restore negative MSBs if needed
//		exponent = exponent << (32 - fifo->exponentSize);
//		exponent >>= (32 - fifo->exponentSize);

		// Substract custom exponent bias
		//exponent += FP32_EXPONENT_BIAS;
//		exponent -= (uint32_t) powf(2,fifo->exponentSize-1)-1;
		exponent -= fifo->exponentBias;
		exponent += FP32_EXPONENT_BIAS;
		exponent &= 0x000000ff;

		// Set exponent position
		exponent <<= (FP32_MANTISSA_SIZE);

		// Get mantissa
		uint32_t mantissa = (*temp >> (fifo->workingSize - fifo->storageSize - (FP32_MANTISSA_SIZE - fifo->mantissaSize))) & 0x007FFFFF;

		// Assemble cFP
		*temp = sign | exponent | mantissa;
	}

//	__asm__("nop");

	return return_value_cFP;
}


void fifoSet_float_cFP(StorageFifo* fifo, float data_in, int index)
{
//	int32_t temp_buffer;

//	uint8_t* data_in_arr = (uint8_t*) &temp_buffer;

//	printf("%x\n", *((uint32_t*) &data_in));

	// Zeroing values too small to be represented
	if (getExponentValue(data_in) + fifo->exponentBias < 0 && data_in != 0)
	{
//		__asm__("nop");
		data_in = 0;
//		uint32_t* temp = (uint32_t*) &data_in;
//		*temp = signbit(data_in) | 0x00000001;
	}

	if (data_in != 0)
	{
		uint32_t* temp = (uint32_t*) &data_in;

		// Get sign
		uint32_t sign = *temp & 0x80000000;

		// Get exponent, on 8 bits or less
		uint32_t exponent = ((uint32_t) getExponentValue(data_in));

		// Apply fifo exponent bias
		exponent += fifo->exponentBias;

		// Crop exponent MSB
		exponent &= ((1 << (fifo->exponentSize))-1);
		// Set exponent position
		exponent <<= (fifo->mantissaSize + (fifo->workingSize - fifo->storageSize));

		// Get mantissa
		uint32_t mantissa = *temp & 0x007FFFFF;
		// Crop mantissa
		mantissa >>= (FP32_MANTISSA_SIZE - fifo->mantissaSize);
		// Set mantissa position
		mantissa <<= (fifo->workingSize - fifo->storageSize);

		// Assemble cFP
		*temp = sign | exponent | mantissa;

//		__asm__("nop");
	}

	// Get int32_t* pointer to store
	int32_t* temp = (int32_t*) &data_in;

	FIFO_SET_FUNC(fifo, *temp, index);
}

void fifoSimulate_float_cFP(int NB_ELEMENT, float* buffer, int storageSize)
{
//	float max_abs_exponent = 0.0f;
	float max_exponent = -INFINITY;
	float min_exponent = INFINITY;
//	float max_abs_exponent_test = 0.0f;
	float error = 0.0f;

	int sign = 0;
//	int sign_test = 0;

	// Get highest value
	for (int i = 0; i < NB_ELEMENT; i++)
	{

		if (buffer[i] != 0)
		{
//			if (abs(getExponentValue(buffer[i])) > max_abs_exponent)
//				max_abs_exponent = abs(getExponentValue(buffer[i]));

			if (getExponentValue(buffer[i]) > max_exponent)
				max_exponent = getExponentValue(buffer[i]);

			if (getExponentValue(buffer[i]) < min_exponent)
				min_exponent = getExponentValue(buffer[i]);

			// Double logical NOT to force non-zero values to 1
//			sign += !(!(signbit(buffer[i])));
		}
	}

//	if (sign == 0)					// All values are positive
//	{
//		sign = 0;
//	}
//	else if (sign == NB_ELEMENT)	// All values are negative
//	{
//		sign = 1;
//	}
//	else
//	{
//		sign = -1;
//	}


//	StorageFifo fifo = fifoCreate_cFP(NB_ELEMENT, storageSize, 32, max_abs_exponent);
	StorageFifo fifo = fifoCreate_cFP_better(NB_ELEMENT, storageSize, 32, min_exponent, max_exponent);


	for (int i = 0; i < NB_ELEMENT; i++)
	{
//		fifoSet_float_cFP(&fifo, buffer[i], i);
		fifoSet_float(&fifo, buffer[i], i);
	}

	for (int i = 0; i < NB_ELEMENT; i++)
	{
//		float temp = fifoGet_float_cFP(&fifo, i);
		float temp = fifoGet_float(&fifo, i);
		error += fabsf(temp - buffer[i]);
		if (fabsf(temp - buffer[i]) != 0)
				__asm__("nop");
		buffer[i] = temp;
	}

	if (error != 0)
	{
		printf("Error = %f\n", error/NB_ELEMENT);
		__asm__("nop");
	}


	fifoDestroy(&fifo);
}


void fifoSet_float_cFP_control(StorageFifo* fifo, float data_in, int index)
{

	int32_t storageIntegerSize = ((int) (floorf(log2f(data_in)) +2 ));

	// Check if data_in can be stored with fifo->
	if (storageIntegerSize > (fifo->storageSize - fifo->storageDecimalSize))
	{
		// if data_in cannot be stored as is

		// determine new fifo->storageDecimalSize but do not change it yet
		int32_t newStorageDecimalSize = fifo->storageSize - storageIntegerSize;

		for (int i = 0; i < fifo->nb_element; i++)
		{
			float temp = fifoGet_float_FxP(fifo, i);
			fifoSet_float_FxP_override(fifo, temp, i, newStorageDecimalSize);

			float check = fifoGet_float_FxP_override(fifo, i, newStorageDecimalSize);
			__asm__("nop");
		}
		fifo->storageDecimalSize = newStorageDecimalSize;
	}

	fifoSet_float_FxP(fifo, data_in, index);
}


// VECTOR generic
float fifoVecGet_float(StorageFifoVec* fifo, int index, int component)
{
	float return_value;

	switch (fifo->fifoType)
	{
		case FxP :
			return_value = fifoVecGet_float_FxP(fifo, index, component);
			break;
		case FxP_control :
			return_value = fifoVecGet_float_FxP(fifo, index, component);
			break;
		case cFP :
			return_value = fifoVecGet_float_cFP(fifo, index, component);
			break;
		case SCALING :
			return_value = fifoVecGet_float_scaling(fifo, index, component);
			break;
		case SCALING_control :
//			return_value = fifoVecGet_float_scaling(fifo, index, component);
			break;
		case NO_TYPE :
		{
			uint32_t temp = fifoVecGet_int32(fifo, index, component);
			return_value = *((float*) &temp);
			break;
		}
		default :
		{
			uint32_t temp = fifoVecGet_int32(fifo, index, component);
			return_value = *((float*) &temp);
		}
	}


	return return_value;
}


float fifoVecSet_float(StorageFifoVec* fifo, float data_in, int index, int component)
{
	switch (fifo->fifoType)
	{
		case FxP :
			fifoVecSet_float_FxP(fifo, data_in, index, component);
			break;
		case FxP_control :
			fifoVecSet_float_FxP_control(fifo, data_in, index, component);
			break;
		case cFP :
			fifoVecSet_float_cFP(fifo, data_in, index, component);
			break;
		case SCALING :
			fifoVecSet_float_scaling(fifo, data_in, index, component);
			break;
		case SCALING_control :
//			fifoVecSet_float_scaling_control(fifo, data_in, index, component);
			break;
		case NO_TYPE :
		{
			uint32_t *temp = (uint32_t*) &data_in;
			fifoVecSet_int32(fifo, *temp, index, component);
			break;
		}
		default :
		{
			uint32_t *temp = (uint32_t*) &data_in;
			fifoVecSet_int32(fifo, *temp, index, component);
		}
	}

	return data_in;
}


uint32_t fifoVecGet_int32(StorageFifoVec* fifo, int index, int component)
{
	uint32_t data_out = 0;

//	int32_t storage_int32_index = index * fifo->storageSize / 32;
//	int32_t storage_int16_index = index * fifo->storageSize / 16;
//	int32_t storage_int8_index = index * fifo->storageSize / 8;

	uint32_t realIndex = index * fifo->nbComponent + component;

//	return FIFO_GET_FUNC(fifo, realIndex);

	uint32_t storage_int32_index = (realIndex * fifo->storageSize) >> 5;
	uint32_t storage_int16_index = (realIndex * fifo->storageSize) >> 4;
	uint32_t storage_int8_index = (realIndex * fifo->storageSize) >> 3;

//	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) % 32);
	uint32_t nb_bit_left = 32 - ((realIndex * fifo->storageSize) & 0x0000001F);

	// nb_bit_right = max(0, blabla) ?
//	int32_t nb_bit_right = max(fifo->storageSize - nb_bit_left, 0);
	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;

	// Left side. There is always a left side

	// Retrieve data from storage (it's inverted)
	uint32_t storage_space_left_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index];

	uint32_t data_to_get_left = BYTE_SWAP32(storage_space_left_inv);

	// Putting MSBs back in place
	data_to_get_left <<= (32 - nb_bit_left);

#if 0
	if (nb_bit_left > fifo->storageSize)
	{
		// Removing bits from the next data
//		data_to_get_left &= ~(0xFFFFFFFF >> (32 - (nb_bit_left - fifo->storageSize)));
		data_to_get_left &= ~(0xFFFFFFFF >> (32 + nb_bit_right));
	}
#else
	// If nb_bit_right >= 0, then no mask is needed. Applying it anyways should do nothing
//	data_to_get_left &= ~(0xFFFFFFFFULL >> (32 + nb_bit_right));
	data_to_get_left &= ~(0xFFFFFFFFULL >> (fifo->storageSize));
#endif

	data_out = data_to_get_left;

	// Right side, if needed
	if (nb_bit_right > 0)
	{
		uint32_t storage_space_right_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index + 1];

		uint32_t data_to_get_right = BYTE_SWAP32(storage_space_right_inv);

		// Putting LSBs back in place
		data_to_get_right >>= (nb_bit_left);

		// Removing bits from the next data
		data_to_get_right &= ~(0xFFFFFFFF >> (fifo->storageSize));

		data_out |= data_to_get_right;
	}

	return data_out;
}

void fifoVecSet_int32(StorageFifoVec* fifo, uint32_t data_in, int index, int component)
{
	// Variable with _inv suffixe are bytes-inverted data, can't perform shift on them, only INVERTED masking
	// Shift operations need to be performed on non-inverted data

//	int32_t storage_int32_index = index * fifo->storageSize / 32;
//	int32_t storage_int16_index = index * fifo->storageSize / 16;
//	int32_t storage_int8_index = index * fifo->storageSize / 8;

	uint32_t realIndex = index * fifo->nbComponent + component;

//	FIFO_SET_FUNC(fifo, data_in, realIndex);
//	return;

	uint32_t storage_int32_index = (realIndex * fifo->storageSize) >> 5;
	uint32_t storage_int16_index = (realIndex * fifo->storageSize) >> 4;
	uint32_t storage_int8_index = (realIndex * fifo->storageSize) >> 3;

//	uint32_t nb_bit_left = 32 - ((index * fifo->storageSize) % 32);
	uint32_t nb_bit_left = 32 - ((realIndex * fifo->storageSize) & 0x0000001F);


//	int32_t nb_bit_right = max(fifo->storageSize - nb_bit_left, 0);
	int32_t nb_bit_right = fifo->storageSize - nb_bit_left;

	// Left side. There is always a left side

	// for a 32bit fetch
	uint32_t storage_space_left_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index];

	// Removing unwanted LSBs
	data_in &= ~(0xFFFFFFFFULL >> (fifo->storageSize));

	// Getting left-side bit to store, then inverting the bytes for storage
	uint32_t data_to_store_left = data_in >> (32 - nb_bit_left);
	uint32_t data_to_store_left_inv = BYTE_SWAP32(data_to_store_left);

	// if all of data_in fit in storage_space_left +1
	#if 0
	if (nb_bit_left > fifo->storageSize)
	{
		storage_space_left_inv = storage_space_left_inv & BYTE_SWAP32(~(0xFFFFFFFF >> (32 - nb_bit_left)) | (0xFFFFFFFF >> (32 - (nb_bit_left - fifo->storageSize))));
	}
	else
	{
		// Prepare new segment
		storage_space_left_inv = storage_space_left_inv & BYTE_SWAP32(~(0xFFFFFFFF >> (32 - nb_bit_left)));
	}
	#else
	storage_space_left_inv = storage_space_left_inv & BYTE_SWAP32(~(0xFFFFFFFF >> (32 - nb_bit_left)) | (0xFFFFFFFFULL >> (32 - (nb_bit_left - fifo->storageSize))));
	#endif

	// Write left segment
	((uint32_t*) fifo->pointerToBuffer)[storage_int32_index] = storage_space_left_inv | data_to_store_left_inv;


	// Right side, if needed
	if (nb_bit_right > 0)
	{
		uint32_t storage_space_right_inv = ((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1];

		// scrapping already stored MSBs
		uint32_t data_to_store_right = (data_in << (nb_bit_left));

		// scrapping unwanted LSBs
		data_to_store_right &= ~((1 << (32 - nb_bit_right)) -1);

		// inverting data for coming insertion
		uint32_t data_to_store_right_inv = BYTE_SWAP32(data_to_store_right);

		// Cleaning destination
		storage_space_right_inv &= BYTE_SWAP32((0x01 << (32 - nb_bit_right)) -1);

		((uint32_t*) fifo->pointerToBuffer)[storage_int32_index +1] = storage_space_right_inv | data_to_store_right_inv;
	}
}

void fifoVecSimulate_float(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent)
{

	float* buffer_array = (float*) buffer;

	float error = 0.0f;

	StorageFifoVec fifo = fifoVecCreate(NB_ELEMENT, storageSize, 32, nbComponent);

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
			fifoVecSet_float(&fifo, buffer_array[j * nbComponent + i], j, i);
		}
//		__asm__("nop");
	}

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
			float temp = fifoVecGet_float(&fifo, j, i);
			error += fabsf(temp - buffer_array[j * nbComponent + i]);
			if (fabsf(temp - buffer_array[j * nbComponent + i]) != 0)
					__asm__("nop");
			buffer_array[j * nbComponent + i] = temp;
		}
	}

	if (error != 0)
		printf("error: %f\n", error);

	fifoVecDestroy(&fifo);
}


// cFP VECTOR

StorageFifoVec fifoVecCreate_cFP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* maxAbsExponent, int nbComponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	StorageFifoVec fifo;

	fifo.storageDecimalSize = NULL;

	fifo.mantissaSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.exponentSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.nbComponent = nbComponent;

	for (int i = 0; i < nbComponent; i++)
	{
		int exponentSize = (int) (floorf(log2f(maxAbsExponent[i])) +2);
		int mantissaSize = storageSize - exponentSize -1;

		if (mantissaSize > FP32_MANTISSA_SIZE)
			mantissaSize = FP32_MANTISSA_SIZE;

		fifo.mantissaSize[i] = mantissaSize;
		fifo.exponentSize[i] = exponentSize;
	}

	fifo.pointerToBuffer = (uint8_t*) malloc(nb_element * storageSize / CHAR_BIT * nbComponent);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;

	return fifo;
}

StorageFifoVec fifoVecCreate_cFP_better(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* minExponent, float* maxExponent, int nbComponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	StorageFifoVec fifo;

	fifo.exponentBias = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.exponentSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.mantissaSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.nbComponent = nbComponent;

	for (int i = 0; i < nbComponent; i++)
	{
		int exponentSize = (int) (floorf(log2f(maxExponent[i] - minExponent[i] +1)+1));
		int32_t exponentBias = -minExponent[i] +1;
		int mantissaSize = storageSize - exponentSize -1;

//		int exponentSize = (int) (floorf(log2f(maxAbsExponent[i])) +2);
//		int mantissaSize = storageSize - exponentSize -1;

		if (mantissaSize > FP32_MANTISSA_SIZE)
			mantissaSize = FP32_MANTISSA_SIZE;

		fifo.exponentBias[i] = exponentBias;
		fifo.mantissaSize[i] = mantissaSize;
		fifo.exponentSize[i] = exponentSize;
	}

	// Division while rounding up
	uint64_t byteToAlloc = (nb_element * nbComponent * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

	//fifo.pointerToBuffer = (uint8_t*) malloc((double) (nb_element * storageSize * nbComponent) / 8.0);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);
//	fifo.pointerToBuffer = (uint8_t*) calloc(byteToAlloc,1);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;

	fifo.fifoType = cFP;

	fifo.storageDecimalSize = NULL;
	fifo.offset = NULL;
	fifo.stepSize = NULL;

	return fifo;
}

StorageFifoVec fifoVecCreate_cFP_better_sign(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* minExponent, float* maxExponent, int* sign, int nbComponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	StorageFifoVec fifo;

	fifo.exponentBias = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.exponentSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.mantissaSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.signBit = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.nbComponent = nbComponent;

	for (int i = 0; i < nbComponent; i++)
	{
		int exponentSize = (int) (floorf(log2f(maxExponent[i] - minExponent[i] +1)+1));
		int32_t exponentBias = -minExponent[i] +1;
		int mantissaSize = storageSize - exponentSize;

//		int exponentSize = (int) (floorf(log2f(maxAbsExponent[i])) +2);
//		int mantissaSize = storageSize - exponentSize -1;

		if (sign[i] == SIGN_UNDEF)
			mantissaSize--;

		if (mantissaSize > FP32_MANTISSA_SIZE)
			mantissaSize = FP32_MANTISSA_SIZE;

		fifo.exponentBias[i] = exponentBias;
		fifo.mantissaSize[i] = mantissaSize;
		fifo.exponentSize[i] = exponentSize;
		fifo.signBit[i] = sign[i];
	}

	// Division while rounding up
	uint64_t byteToAlloc = (nb_element * nbComponent * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

	//fifo.pointerToBuffer = (uint8_t*) malloc((double) (nb_element * storageSize * nbComponent) / 8.0);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);
//	fifo.pointerToBuffer = (uint8_t*) calloc(byteToAlloc,1);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;

	fifo.fifoType = cFP;

	fifo.storageDecimalSize = NULL;
	fifo.offset = NULL;
	fifo.stepSize = NULL;

	return fifo;
}

float fifoVecGet_float_cFP(StorageFifoVec* fifo, int index, int component)
{
	float return_value_cFP = 0.0f;

	// fifoGet_int32 return data that is right-padded with 0s
	int32_t* return_value = (int32_t*) &return_value_cFP;

	*return_value = fifoVecGet_int32(fifo, index, component);

	if (return_value_cFP != 0)
	{
		// Get pointer
		uint32_t* temp = (uint32_t*) &return_value_cFP;

		// Get sign
		uint32_t sign = *temp & 0x80000000;

		// Get exponent, on 8 bits or less
		int32_t exponent = (*temp >> (fifo->mantissaSize[component] + (fifo->workingSize - fifo->storageSize)));

		// Isolate exponent bits
		exponent &= ((1<<(fifo->exponentSize[component]))-1);

		// Restore negative MSBs if needed
//		exponent = exponent << (32 - fifo->exponentSize);
//		exponent >>= (32 - fifo->exponentSize);

		// Substract custom exponent bias
		//exponent += FP32_EXPONENT_BIAS;
//		exponent -= (uint32_t) powf(2,fifo->exponentSize-1)-1;
		exponent -= fifo->exponentBias[component];
		exponent += FP32_EXPONENT_BIAS;
		exponent &= 0x000000ff;

		// Set exponent position
		exponent <<= (FP32_MANTISSA_SIZE);

		// Get mantissa
		uint32_t mantissa = (*temp >> (fifo->workingSize - fifo->storageSize - (FP32_MANTISSA_SIZE - fifo->mantissaSize[component]))) & 0x007FFFFF;

		// Assemble cFP
		*temp = sign | exponent | mantissa;
	}

//	__asm__("nop");

	return return_value_cFP;
}

void fifoVecSet_float_cFP(StorageFifoVec* fifo, float data_in, int index, int component)
{
	//	int32_t temp_buffer;

//	uint8_t* data_in_arr = (uint8_t*) &temp_buffer;

//	printf("%x\n", *((uint32_t*) &data_in));

	// Zeroing values too small to be represented
	if (getExponentValue(data_in) + fifo->exponentBias[component] < 0 && data_in != 0)
	{
//		__asm__("nop");
		data_in = 0;
//		uint32_t* temp = (uint32_t*) &data_in;
//		*temp = signbit(data_in) | 0x00000001;
	}

	if (data_in != 0)
	{
		uint32_t* temp = (uint32_t*) &data_in;

		// Get sign
		uint32_t sign = *temp & 0x80000000;

		// Get exponent, on 8 bits or less
		uint32_t exponent = ((uint32_t) getExponentValue(data_in));

		// Apply fifo exponent bias
		exponent += fifo->exponentBias[component];

		// Crop exponent MSB
		exponent &= ((1 << (fifo->exponentSize[component]))-1);
		// Set exponent position
		exponent <<= (fifo->mantissaSize[component] + (fifo->workingSize - fifo->storageSize));

		// Get mantissa
		uint32_t mantissa = *temp & 0x007FFFFF;
		// Crop mantissa
		mantissa >>= (FP32_MANTISSA_SIZE - fifo->mantissaSize[component]);
		// Set mantissa position
		mantissa <<= (fifo->workingSize - fifo->storageSize);

		// Assemble cFP
		*temp = sign | exponent | mantissa;

//		__asm__("nop");
	}

	// Get int32_t* pointer to store
	int32_t* temp = (int32_t*) &data_in;

	fifoVecSet_int32(fifo, *temp, index, component);
}

void fifoVecSimulate_float_cFP(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent)
{
	float* max_abs_exponent = (float*) calloc(nbComponent, sizeof(float));
	float* min_exponent = (float*) malloc(nbComponent * sizeof(float));
	float* max_exponent = (float*) malloc(nbComponent * sizeof(float));

	int* sign = (int*) calloc(nbComponent, sizeof(int));

	float* buffer_array = (float*) buffer;

	float error = 0.0f;

	//Get highest values
	for (int i = 0; i < nbComponent; i++)
	{
		min_exponent[i] = FLT_MAX;
		max_exponent[i] = -FLT_MAX;
		for (int j = 0; j < NB_ELEMENT; j++)
		{
			if (buffer_array[j * nbComponent + i] != 0)
			{
				if (abs(getExponentValue(buffer_array[j * nbComponent + i])) > max_abs_exponent[i])
					max_abs_exponent[i] = abs(getExponentValue(buffer_array[j * nbComponent + i]));

				if (getExponentValue(buffer_array[j * nbComponent + i]) > max_exponent[i])
					max_exponent[i] = getExponentValue(buffer_array[j * nbComponent + i]);

				if (getExponentValue(buffer_array[j * nbComponent + i]) < min_exponent[i])
					min_exponent[i] = getExponentValue(buffer_array[j * nbComponent + i]);

				// Double logical NOT to force non-zero values to 1
				sign[i] += !(!(signbit(buffer_array[j * nbComponent + i])));
			}
		}
	}

	for (int i = 0; i < nbComponent; i++)
	{
//		uint32_t sign_test = 0;
		if (sign[i] == 0)					// All values are positive
		{
			sign[i] = SIGN_POS;
		}
		else if (sign[i] == NB_ELEMENT)	// All values are negative
		{
			sign[i] = SIGN_NEG;
		}
		else
		{
			sign[i] = SIGN_UNDEF;
		}
	}

	//StorageFifoVec fifo = fifoVecCreate_cFP(NB_ELEMENT, storageSize, 32, max_abs_exponent, nbComponent);
	StorageFifoVec fifo = fifoVecCreate_cFP_better(NB_ELEMENT, storageSize, 32, min_exponent, max_exponent, nbComponent);

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
			fifoVecSet_float_cFP(&fifo, buffer_array[j * nbComponent + i], j, i);
		}
		__asm__("nop");
	}

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
			float temp = fifoVecGet_float_cFP(&fifo, j, i);

			error += fabsf(temp - buffer_array[j * nbComponent + i]);
			if (fabsf(temp - buffer_array[j * nbComponent + i]) != 0)
					__asm__("nop");
			buffer_array[j * nbComponent + i] = temp;

			//buffer_array[j * nbComponent + i] = fifoVecGet_float_cFP(&fifo, j, i);
		}
	}

	if (error != 0)
		printf("error: %f\n", error);

	fifoVecDestroy(&fifo);
//	free(sign);
	free(min_exponent);
	free(max_exponent);
	free(max_abs_exponent);
}


double fifoGet_double_FxP(StorageFifo* fifo, int index)
{
	double return_value = 0.0f;
	int64_t return_value_FxP = 0;

	uint8_t* data_out_arr = (uint8_t*) &return_value_FxP;

	uint32_t storageByteIndex = index * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = index * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;
	uint32_t remainingBitToRetrieve = fifo->storageSize;

	uint32_t i = 0;

	while(remainingBitToRetrieve > CHAR_BIT)
	{
		data_out_arr[workingByteIndex - (i+1)] =  (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | (fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment);
		i++;
		remainingBitToRetrieve -= CHAR_BIT;
	}

	// cas spécifique dernier/seul octet
	if (((index+1)*fifo->storageSize  & 0x07) >= (remainingBitToRetrieve & 0x07))
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1);
	}
	else
	{
		data_out_arr[workingByteIndex - (i+1)] = ((fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | ((fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1)));
	}

	return_value_FxP = return_value_FxP >> (fifo->workingSize - fifo->storageSize);

	return_value = (double) return_value_FxP;

	#if 0

	if (return_value != 0)
	{
		uint64_t* temp = (uint64_t*) &return_value;

		uint8_t exponent = (*temp >> 23) & 0xff;
		exponent -= fifo->storageDecimalSize;

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << 23);
	}

	#else

	return_value = return_value / pow(2, fifo->storageDecimalSize);

	#endif

	__asm__("nop");

	return return_value;
}


void fifoSet_double_FxP(StorageFifo* fifo, double data_in, int index)
{
	int32_t temp_buffer;

	uint8_t* data_in_arr = (uint8_t*) &temp_buffer;

//	printf("%x\n", *((uint32_t*) &data_in));

#if 0

	if (data_in != 0)
	{
		uint32_t* temp = (uint32_t*) &data_in;

		uint8_t exponent = (*temp >> 23) & 0xff;

		exponent += fifo->storageDecimalSize;

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << 23);
	}

#else

	data_in = data_in * pow(2, fifo->storageDecimalSize);

#endif


//	printf("%x\n", *((uint32_t*) &data_in));

	temp_buffer = (int64_t) data_in;

	temp_buffer = temp_buffer << (fifo->workingSize - fifo->storageSize);

	uint32_t storageByteIndex = index * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = index * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;

	int i = 0;

	uint32_t remainingBitToStore = fifo->storageSize;

	uint32_t widowBits = min(remainingBitToStore, rightMisalignment & 0x07);

	// If widow bits
	if (widowBits != 0)
	{
		fifo->pointerToBuffer[storageByteIndex] = ((data_in_arr[workingByteIndex - 1] >> leftMisalignment) & ((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment))
												  | (fifo->pointerToBuffer[storageByteIndex] & (uint8_t) ~((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment));
		i++;
		remainingBitToStore -= widowBits;
	}

	int j=1;

	while(remainingBitToStore > CHAR_BIT)
	{
		fifo->pointerToBuffer[storageByteIndex + i++] = (data_in_arr[workingByteIndex - j] << (rightMisalignment & 0x07))
														| (data_in_arr[workingByteIndex - j - 1] >> (CHAR_BIT - (rightMisalignment & 0x07)));
		j++;
		remainingBitToStore -= CHAR_BIT;
	}

	// cas spécifique dernier octet
	if (remainingBitToStore != 0)
	{
		if (widowBits + remainingBitToStore <= CHAR_BIT)			// If orphan bits are on the same byte
		{
			fifo->pointerToBuffer[storageByteIndex + i] = ((data_in_arr[workingByteIndex - j] << (rightMisalignment  & 0x07)) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
		else												// If orphan bits are on 2 bytes
		{
			fifo->pointerToBuffer[storageByteIndex + i] = (((data_in_arr[workingByteIndex - j] << rightMisalignment) | data_in_arr[workingByteIndex - j - 1] >> leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
	}
}


StorageFifoVec fifoVecCreate(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, uint8_t nbComponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	StorageFifoVec fifo;

	fifo.nbComponent = nbComponent;

	fifo.pointerToBuffer = (uint8_t*) malloc((nb_element * storageSize / CHAR_BIT) * nbComponent);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;

	fifo.storageDecimalSize = NULL;

	fifo.exponentBias = NULL;
	fifo.exponentSize = NULL;
	fifo.mantissaSize = NULL;

	fifo.offset = NULL;
	fifo.stepSize = NULL;

	return fifo;
}



void fifoVecDestroy(StorageFifoVec* fifo)
{
	if (fifo->pointerToBuffer != NULL)
	{
		free(fifo->pointerToBuffer);
		fifo->pointerToBuffer = NULL;
	}

	// fxp
	if (fifo->storageDecimalSize != NULL)
	{
		free(fifo->storageDecimalSize);
		fifo->storageDecimalSize = NULL;
	}

	// cfp
	if (fifo->exponentBias != NULL)
	{
		free(fifo->exponentBias);
		fifo->exponentBias = NULL;
	}
	if (fifo->exponentSize != NULL)
	{
		free(fifo->exponentSize);
		fifo->exponentSize = NULL;
	}
	if (fifo->mantissaSize != NULL)
	{
		free(fifo->mantissaSize);
		fifo->mantissaSize = NULL;
	}

	// scaling
	if (fifo->offset == NULL)
	{
		free(fifo->offset);
		fifo->offset = NULL;
	}
	if (fifo->stepSize == NULL)
	{
		free(fifo->stepSize);
		fifo->stepSize = NULL;
	}
}

float fifoVecGet_float_FxP(StorageFifoVec* fifo, int index, int component)
{
	float return_value = 0.0f;
	#if 0
	int32_t return_value_FxP = 0;

	uint8_t* data_out_arr = (uint8_t*) &return_value_FxP;

	uint32_t realIndex = index * fifo->nbComponent + component;

	uint32_t storageByteIndex = realIndex * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = realIndex * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;
	uint32_t remainingBitToRetrieve = fifo->storageSize;

	uint32_t i = 0;

	while(remainingBitToRetrieve > CHAR_BIT)
	{
		data_out_arr[workingByteIndex - (i+1)] =  (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | (fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment);
		i++;
		remainingBitToRetrieve -= CHAR_BIT;
	}

	// cas spécifique dernier/seul octet
	if (((realIndex+1)*fifo->storageSize  & 0x07) >= (remainingBitToRetrieve & 0x07))
	{
		data_out_arr[workingByteIndex - (i+1)] = (fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1);
	}
	else
	{
		data_out_arr[workingByteIndex - (i+1)] = ((fifo->pointerToBuffer[storageByteIndex + i] << leftMisalignment)
												  | ((fifo->pointerToBuffer[storageByteIndex + 1 + i] >> rightMisalignment) & ~((1 << (CHAR_BIT - remainingBitToRetrieve)) - 1)));
	}

	return_value_FxP = return_value_FxP >> (fifo->workingSize - fifo->storageSize);

	return_value = (float) return_value_FxP;
	#else

	int32_t temp = fifoVecGet_int32(fifo, index, component);

	temp >>= (fifo->workingSize - fifo->storageSize);

	return_value = (float) temp;

	#endif

	if (return_value != 0)
	{
		uint32_t* temp = (uint32_t*) &return_value;

		uint8_t exponent = (*temp >> FP32_MANTISSA_SIZE) & 0xff;
		exponent -= fifo->storageDecimalSize[component];

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << FP32_MANTISSA_SIZE);
	}

	__asm__("nop");

	return return_value;

}


void fifoVecSet_float_FxP(StorageFifoVec* fifo, float data_in, int index, int component)
{

//	printf("%x\n", *((uint32_t*) &data_in));

	if (data_in != 0)
	{
		uint32_t* temp = (uint32_t*) &data_in;

		uint8_t exponent = (*temp >> FP32_MANTISSA_SIZE) & 0xff;

		exponent += fifo->storageDecimalSize[component];

		*temp = (*temp & 0x807fffff) | (((uint32_t) exponent) << FP32_MANTISSA_SIZE);
	}

	#if 1

	// Casting data_in to int32 to complete FxP conversion
	int32_t temp = (int32_t) data_in;

	// Setting temp as right-padded to work with fifoSet_int32
	temp <<= (fifo->workingSize - fifo->storageSize);

	fifoVecSet_int32(fifo, temp, index, component);

	#else

	int32_t temp_buffer;

	uint8_t* data_in_arr = (uint8_t*) &temp_buffer;

//	printf("%x\n", *((uint32_t*) &data_in));

	temp_buffer = (int32_t) data_in;

	temp_buffer = temp_buffer << (fifo->workingSize - fifo->storageSize);

	uint32_t realIndex = index * fifo->nbComponent + component;

	uint32_t storageByteIndex = realIndex * fifo->storageSize >> 3;
	uint32_t workingByteIndex = fifo->workingSize >> 3;
	uint8_t leftMisalignment = realIndex * fifo->storageSize & 0x07;
	uint8_t rightMisalignment = CHAR_BIT - leftMisalignment;

	int i = 0;

	uint32_t remainingBitToStore = fifo->storageSize;

	uint32_t widowBits = min(remainingBitToStore, rightMisalignment & 0x07);

	// If widow bits
	if (widowBits != 0)
	{
		fifo->pointerToBuffer[storageByteIndex] = ((data_in_arr[workingByteIndex - 1] >> leftMisalignment) & ((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment))
												  | (fifo->pointerToBuffer[storageByteIndex] & (uint8_t) ~((uint8_t) ~((1 << (CHAR_BIT - widowBits)) - 1) >> leftMisalignment));
		i++;
		remainingBitToStore -= widowBits;
	}

	int j=1;

	while(remainingBitToStore > CHAR_BIT)
	{
		fifo->pointerToBuffer[storageByteIndex + i++] = (data_in_arr[workingByteIndex - j] << (rightMisalignment & 0x07))
														| (data_in_arr[workingByteIndex - j - 1] >> (CHAR_BIT - (rightMisalignment & 0x07)));
		j++;
		remainingBitToStore -= CHAR_BIT;
	}

	// cas spécifique dernier octet
	if (remainingBitToStore != 0)
	{
		if (widowBits + remainingBitToStore <= CHAR_BIT)			// If orphan bits are on the same byte
		{
			fifo->pointerToBuffer[storageByteIndex + i] = ((data_in_arr[workingByteIndex - j] << (rightMisalignment  & 0x07)) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
		else												// If orphan bits are on 2 bytes
		{
			fifo->pointerToBuffer[storageByteIndex + i] = (((data_in_arr[workingByteIndex - j] << rightMisalignment) | data_in_arr[workingByteIndex - j - 1] >> leftMisalignment) & ~((1 << (CHAR_BIT - remainingBitToStore)) - 1))
														  | (fifo->pointerToBuffer[storageByteIndex + i] & ((1 << (CHAR_BIT - remainingBitToStore)) - 1));
		}
	}

	#endif

}

void fifoVecSet_float_FxP_control(StorageFifoVec* fifo, float data_in, int index, int component)
{
	#warning fifoVecSet_float_FxP_control not implemented yet
}


StorageFifoVec fifoVecCreate_FxP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* maxValue, int nbComponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

	StorageFifoVec fifo;

	fifo.storageDecimalSize = (int32_t*) malloc(nbComponent * sizeof(int32_t));
	fifo.nbComponent = nbComponent;

	for (int i = 0; i < nbComponent; i++)
	{
		int storageIntegerSize = (int) (floorf(log2f(maxValue[i])) +2);
		int storageDecimalSize = storageSize - storageIntegerSize;

		fifo.storageDecimalSize[i] = storageDecimalSize;
	}

	uint64_t byteToAlloc = (nb_element * nbComponent * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

	//fifo.pointerToBuffer = (uint8_t*) malloc((double) (nb_element * storageSize * nbComponent) / 8.0);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;

	fifo.exponentSize = NULL;
	fifo.exponentBias = NULL;
	fifo.mantissaSize = NULL;
	fifo.signBit = NULL;

	fifo.offset = NULL;
	fifo.stepSize = NULL;

	fifo.fifoType = FxP;

	return fifo;
}



void fifoVecSimulate_float_FxP(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent)
{
	float* max_value = (float*) calloc(nbComponent, sizeof(float));

	float* buffer_array = (float*) buffer;

	float error = 0;

	//Get highest values
	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
			if (fabsf(buffer_array[j * nbComponent + i]) > max_value[i])
			{
//				max_value[i] = fabsf(buffer_array[j * nbComponent + i]);
				max_value[i] = buffer_array[j * nbComponent + i];
			}
		}
	}

	StorageFifoVec fifo = fifoVecCreate_FxP(NB_ELEMENT, storageSize, 32, max_value, nbComponent);

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
//			fifoVecSet_float_FxP(&fifo, buffer_array[j * nbComponent + i], j, i);
			fifoVecSet_float(&fifo, buffer_array[j * nbComponent + i], j, i);
		}
		__asm__("nop");
	}

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
//			float temp = fifoVecGet_float_FxP(&fifo, j, i);
			float temp = fifoVecGet_float(&fifo, j, i);
//			printf("%f\t", buffer_array[j * nbComponent + i]);
			if (temp != buffer_array[j * nbComponent + i])
				error += fabsf(temp - buffer_array[j * nbComponent + i]);

			buffer_array[j * nbComponent + i] = temp;
//			printf("%f\n", buffer_array[j * nbComponent + i]);
		}
	}

	printf("Avg error = %f\n", error/NB_ELEMENT);

	fifoVecDestroy(&fifo);
	free(max_value);
}

// scaling VECTOR
StorageFifoVec fifoVecCreate_scaling(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* minValue, float* maxValue, int nbComponent)
{
	if (storageSize > workingSize)
	{
		printf("Can't create fifo in %s:%d\n", __func__, __LINE__);
		exit(EXIT_FAILURE);
	}

//	float amplitude = maxValue - minValue;

	StorageFifoVec fifo;

	// Division while rounding up
	uint64_t byteToAlloc = (nb_element * nbComponent * storageSize + (CHAR_BIT-1)) / CHAR_BIT;

	//fifo.pointerToBuffer = (uint8_t*) malloc((double) (nb_element * storageSize * nbComponent) / 8.0);
	fifo.pointerToBuffer = (uint8_t*) malloc(byteToAlloc);

	fifo.storageSize = storageSize;
	fifo.workingSize = workingSize;
	fifo.nb_element = nb_element;
	fifo.nbComponent = nbComponent;

	fifo.fifoType = SCALING;
//	fifo.fifoType = SCALING_control;

	fifo.offset = (float*) malloc(nbComponent * sizeof(float));
	fifo.stepSize = (float*) malloc(nbComponent * sizeof(float));

	for (int i = 0; i < nbComponent; i++)
	{
//		fifo.offset[i] = -minValue[i];
		fifo.stepSize[i] = (maxValue[i] - minValue[i]) * SCALING_AMPLITUDE_CORRECTION / (powf(2, storageSize));
		fifo.offset[i] = (int32_t) (-minValue[i]/fifo.stepSize[i]);
	}

//	fifo.originalStepSize = fifo.stepSize;
//	fifo.originalOffset = fifo.offset;

//	fifo.minValue = minValue;
//	fifo.maxValue = maxValue;

//	fifo.currentMaxValue = -INFINITY;
//	fifo.currentMinValue = INFINITY;

//	fifo.indexOffset = 0;
//	fifo.nbInstance = 1;

	// fxp
	fifo.storageDecimalSize = NULL;

	// cfp
	fifo.exponentBias = NULL;
	fifo.exponentSize = NULL;
	fifo.mantissaSize = NULL;


	return fifo;
}


float fifoVecGet_float_scaling(StorageFifoVec* fifo, int index, int component)
{
	float return_value = 0.0f;

	// fifoGet_int32 return data that is right-padded with 0s
	return_value = (float) ((uint32_t) fifoVecGet_int32(fifo, index, component) >> (fifo->workingSize - fifo->storageSize));

//	return_value = return_value * fifo->stepSize[component] - fifo->offset[component];

	return_value = (return_value - fifo->offset[component]) * fifo->stepSize[component];

	__asm__("nop");

	return return_value;
}


void fifoVecSet_float_scaling(StorageFifoVec* fifo, float data_in, int index, int component)
{
	// Applying scaling
//	data_in = (data_in + fifo->offset[component]) / fifo->stepSize[component];

	data_in = data_in / fifo->stepSize[component] + fifo->offset[component];

	// Casting data_in to int32 to complete FxP conversion
	uint32_t temp = (uint32_t) data_in;

	// Dirty patch, need to actually fix the stepsize
	if (temp >= ((uint64_t) 0x01 << fifo->storageSize))
	{
		__asm__("nop");
//		temp = (0x01 << fifo->storageSize) -1;
		printf("Overflow detected, please raise max value\n");
	}

	// Setting temp as right-padded to work with fifoSet_int32
	temp <<= (fifo->workingSize - fifo->storageSize);

	fifoVecSet_int32(fifo, temp, index, component);
}

void fifoVecSimulate_float_scaling(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent)
{
	float* min_value = (float*) malloc(nbComponent * sizeof(float));
	float* max_value = (float*) malloc(nbComponent * sizeof(float));

	float* buffer_array = (float*) buffer;

	float error = 0;

	//Get highest values
	for (int i = 0; i < nbComponent; i++)
	{
		max_value[i] = -INFINITY;
		min_value[i] = INFINITY;

		for (int j = 0; j < NB_ELEMENT; j++)
		{
			if (buffer_array[j * nbComponent + i] > max_value[i])
				max_value[i] = buffer_array[j * nbComponent + i];

			if (buffer_array[j * nbComponent + i] < min_value[i])
				min_value[i] = buffer_array[j * nbComponent + i];
		}
	}

	StorageFifoVec fifo = fifoVecCreate_scaling(NB_ELEMENT, storageSize, 32, min_value, max_value, nbComponent);

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
//			fifoVecSet_float_FxP(&fifo, buffer_array[j * nbComponent + i], j, i);
			fifoVecSet_float(&fifo, buffer_array[j * nbComponent + i], j, i);
		}
		__asm__("nop");
	}

	for (int i = 0; i < nbComponent; i++)
	{
		for (int j = 0; j < NB_ELEMENT; j++)
		{
//			float temp = fifoVecGet_float_FxP(&fifo, j, i);
			float temp = fifoVecGet_float(&fifo, j, i);
//			printf("%f\t", buffer_array[j * nbComponent + i]);
			if (temp != buffer_array[j * nbComponent + i])
				error += fabsf(temp - buffer_array[j * nbComponent + i]);

			buffer_array[j * nbComponent + i] = temp;
//			printf("%f\n", buffer_array[j * nbComponent + i]);
		}
	}

	printf("Avg error = %f\n", error/NB_ELEMENT);

	fifoVecDestroy(&fifo);
	free(max_value);
	free(min_value);
}

#pragma clang diagnostic pop