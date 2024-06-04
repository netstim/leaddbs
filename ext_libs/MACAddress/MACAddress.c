/* Windows version of C code to get MAC address.
 * To compile in Matlab/Octave:
 *  mex -liphlpapi MACAddress.c
 * will create MACAddress.mex* for the corresponding language. 
 * The mex file, if exists, will take precedence over MACAddress.m.
 *  mac1 = MACAddress(); % return the first MAC address as char;
 *  macAll = MACAddress(1); % return all MAC address as cellstr.
 *  [macAll, st] = MACAddress(1); % also return more info in st.
 *  
 * 170617 Wrote it by Xiangrui.Li at gmail.com 
 * 171029 Remove NICs limit. Make it work for Octave
 * 180627 Rename to the same name as m file (mex takes precedence).
 *        Add 2nd output, switch to GetAdaptersAddresses for IPv6
 */

#include <mex.h>
#include <stdio.h>
#include <winsock2.h>
#include <IPHlpApi.h>
//#pragma comment(lib, "IPHLPAPI.lib")

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  
    char ch18[18], ch[128];
    UINT8 *p, j;
    mwSize i = 0, nNIC = 0;
    ULONG outBufLen = 15000;
    PIP_ADAPTER_ADDRESSES pAddresses = NULL, pCur = NULL;
    PIP_ADAPTER_UNICAST_ADDRESS pUnicast = NULL;
    BOOL allMAC = FALSE;    
    if (nrhs>0) allMAC = (BOOL) mxGetScalar(prhs[0]);
 
    for (j = 0; j<3; j++) {
        pAddresses = (IP_ADAPTER_ADDRESSES *) mxMalloc(outBufLen);
        if (GetAdaptersAddresses(AF_UNSPEC, GAA_FLAG_INCLUDE_PREFIX, NULL, 
                pAddresses, &outBufLen) != ERROR_BUFFER_OVERFLOW) break;
        mxFree(pAddresses);
        pAddresses = NULL;
    }
    if (j>2) { mxFree(pAddresses); pAddresses = NULL;} 
    
    if (allMAC) {
        pCur = pAddresses;
        while (pCur) { // get nNIC only
            if (pCur->PhysicalAddressLength == 6) nNIC++;
            pCur = pCur->Next;
        }
        if (nNIC<1) nNIC = 1;
        plhs[0] = mxCreateCellMatrix(1, nNIC);
    }else nNIC = 1;
            
    if (nlhs>1) {
        const char *field_names[] = {"FriendlyName", "Description", 
                        "MAC_address", "IPv4_address", "IPv6_address"};
        plhs[1] = mxCreateStructMatrix(1, nNIC, 5, field_names);
    }

    for (i=0, pCur=pAddresses; i<nNIC && pCur!=NULL; i++) {
        while (pCur->PhysicalAddressLength != 6) pCur = pCur->Next;
        p = (void*) pCur->PhysicalAddress;
        sprintf(ch18, "%02X", p[0]); // unix %02x
        for (j=1; j<6; j++) sprintf(ch18+strlen(ch18), "-%02X", p[j]);
        if (allMAC) mxSetCell(plhs[0], i, mxCreateString(ch18));
        else plhs[0] = mxCreateString(ch18);
        
        if (nlhs>1) {
            sprintf(ch, "%wS", pCur->FriendlyName);
            mxSetFieldByNumber(plhs[1], i, 0, mxCreateString(ch));
            sprintf(ch, "%wS", pCur->Description);
            mxSetFieldByNumber(plhs[1], i, 1, mxCreateString(ch));
            mxSetFieldByNumber(plhs[1], i, 2, mxCreateString(ch18));
            
            pUnicast = pCur->FirstUnicastAddress;
            while (pUnicast) {
                p = (void*)pUnicast->Address.lpSockaddr; // use hard-coded index
                if (pUnicast->Address.lpSockaddr->sa_family == AF_INET) {
                    sprintf(ch, "%i.%i.%i.%i", p[4], p[5], p[6], p[7]);
                    mxSetFieldByNumber(plhs[1], i, 3, mxCreateString(ch));
                }else if (pUnicast->Address.lpSockaddr->sa_family == AF_INET6) {
                    strcpy(ch, "fe80:"); // fe80:0:0:0: 64 bits for link-local 
                    for (j=16; j<24; j++) { // rest 64 bits
                        if (j%2==0 && p[j] ==0) continue;
                        else if (j%2==1 && p[j-1]!=0) sprintf(ch+strlen(ch), "%02x", p[j]);
                        else sprintf(ch+strlen(ch), ":%x",  p[j]);
                    }
                    mxSetFieldByNumber(plhs[1], i, 4, mxCreateString(ch));
                }
                pUnicast = pUnicast->Next;
            }
        }
        
        pCur = pCur->Next;
    }
    
    if (pAddresses) mxFree(pAddresses);
    else {
        mexWarnMsgIdAndTxt("MACAddress:RandomMAC", "Returned MAC are random numbers");
        sprintf(ch18, "%02X", (UINT8)rand());
        for (j=0; j<5; j++) sprintf(ch18+strlen(ch18), "-%02X", (UINT8)rand());
        if (allMAC) mxSetCell(plhs[0], 0, mxCreateString(ch18));
        else plhs[0] = mxCreateString(ch18);
        if (nlhs>1) {
            mxSetFieldByNumber(plhs[1], 0, 0, mxCreateString("Failed to find network adapter"));
            mxSetFieldByNumber(plhs[1], 0, 1, mxCreateString(ch18));
        }
    }
}
