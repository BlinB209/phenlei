LIB_EXPORT inline FieldProxy * Gradient::GetGradX() const
{
	return dqdx;
}

LIB_EXPORT inline FieldProxy * Gradient::GetGradY() const
{
	return dqdy;
}

LIB_EXPORT inline FieldProxy * Gradient::GetGradZ() const
{
	return dqdz;
}

LIB_EXPORT inline FieldProxy * GradientCellCenter::GetGradX() const
{
	return dqdx;
}

LIB_EXPORT inline FieldProxy * GradientCellCenter::GetGradY() const
{
	return dqdy;
}

LIB_EXPORT inline FieldProxy * GradientCellCenter::GetGradZ() const
{
	return dqdz;
}
