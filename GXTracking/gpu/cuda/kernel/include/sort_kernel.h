#ifndef SORT_KERNEL_H
#define SORT_KERNEL_H

void 
count_by_process_gpu(G4int nTracks, GXTrack *itracks, 
		     G4int *nbrem, G4int *nioni, 
		     int blocksPerGrid, int threadsPerBlock,
                     cudaStream_t stream);

void 
count_by_process_cpu(G4int nTracks, GXTrack *itracks,
		     G4int *nbrem, G4int *nioni);

void 
sort_by_process_gpu(G4int nTracks, GXTrack *itracks, GXTrack *stracks, 
		    G4int *nbrem, G4int *stackSize_brem,
		    G4int *nioni, G4int *stackSize_ioni,
		    int blocksPerGrid, int threadsPerBlock,
                    cudaStream_t stream);

void 
sort_by_process_cpu(G4int nTracks, GXTrack *itracks, GXTrack *stracks,
		    G4int *nbrem, G4int *stackSize_brem,
		    G4int *nioni, G4int *stackSize_ioni);

#endif
