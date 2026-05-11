<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text"/>
  
  <xsl:template match="/">
    <xsl:apply-templates select="//Abstract"/>
    <xsl:apply-templates select="//AbstractText"/>
    <xsl:apply-templates select="//MeshHeading"/>
    <xsl:apply-templates select="//Introduction"/>
    <xsl:apply-templates select="//Methods"/>
    <xsl:apply-templates select="//Results"/>
    <xsl:apply-templates select="//Conclusion"/>
  </xsl:template>

  <xsl:template match="Abstract | AbstractText | MeshHeading | Introduction | Methods | Results | Conclusion">
    <xsl:value-of select="."/>
    <xsl:text>&#10;</xsl:text>
  </xsl:template>
</xsl:stylesheet>
