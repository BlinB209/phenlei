LIB_EXPORT inline void Limiter::SetLimiterTypeID(int uns_limiter) 
{ 
	this->limiterTypeID = uns_limiter; 
}

LIB_EXPORT inline void Limiter::SetLimiterModel(int limit_mode)
{ 
	this->limitModel = limit_mode; 
}

LIB_EXPORT inline void Limiter::SetLimiterVector(int limit_vec) 
{ 
	this->limitVector = limit_vec;
}

LIB_EXPORT inline void Limiter::SetVencatLimiterType(int ivencat) 
{
	this->vencatLimiterType = ivencat; 
}