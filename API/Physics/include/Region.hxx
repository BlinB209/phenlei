inline Zone * Region::GetZone(int iZone)
{
    return (*zoneContainer)[iZone];
}

inline bool Region::ZoneBelongToCurrentProcessor( int globalZoneIndex )
{
    if ( this->GetZone(globalZoneIndex) ) 
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline int Region::GetProcessLocalZoneIndexToGlobalZoneIndex( int localZoneIndex )
{ 
    //return processLocalZoneIndexToGlobalZoneIndex[ localZoneIndex ];
    return PHMPI::GetLocalZoneIDToGlobalZoneID(localZoneIndex);
}